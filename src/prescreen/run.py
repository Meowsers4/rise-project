"""Stage 2 -- cheap empirical prescreen (FoldX + Rosetta cartesian_ddg).

Runs a fast empirical ΔΔG on every variant as (a) a first estimate, (b) a
prioritization ranking, and (c) an independent cross-check against FEP later
(README §4 Stage 2). Both methods report **positive = destabilizing**, matching
``validation.ddg_sign`` so signs line up with the experimental controls and FEP.

FoldX and Rosetta are licensed, SCC-module tools, so the binary invocation is
isolated behind a small backend adapter. The output *parsing*, *aggregation*, and
*controls cross-check* are pure functions, unit-tested with mock tool output; a
:class:`MockBackend` lets the whole path run locally without the binaries.

All parameters come from ``config/pipeline.yaml`` (CLAUDE.md #4).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from src.prep.build import load_config, load_variant_record
from src.seeds import stable_seed

ROOT = Path(__file__).resolve().parents[2]


# --------------------------------------------------------------------------- #
# Output parsers (pure)                                                         #
# --------------------------------------------------------------------------- #
def parse_foldx_ddg(dif_fxout_text: str) -> float:
    """Total ΔΔG (kcal/mol, + = destabilizing) from a FoldX ``Dif_*.fxout`` file.

    The data rows start with the mutated PDB name; the total ΔΔG is the next column.
    """
    for line in dif_fxout_text.splitlines():
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) >= 2 and parts[0].endswith(".pdb"):
            try:
                return float(parts[1])
            except ValueError:
                continue
    raise ValueError("no FoldX data row (line starting with '<pdb>.pdb') found")


def parse_rosetta_cartesian(ddg_text: str, reu_per_kcal: float) -> float:
    """ΔΔG (kcal/mol) from Rosetta ``cartesian_ddg`` output.

    ΔΔG = mean(MUT total REU) - mean(WT total REU), converted REU -> kcal/mol.
    Lines look like ``COMPLEX:  Round<i>:  <WT|MUT>_<...>  <total>  <terms...>``.
    """
    wt, mut = [], []
    for line in ddg_text.splitlines():
        parts = line.split()
        if len(parts) < 4 or not parts[0].startswith("COMPLEX"):
            continue
        label = parts[2].upper()
        try:
            total = float(parts[3])
        except ValueError:
            continue
        if label.startswith("WT"):
            wt.append(total)
        elif label.startswith("MUT"):
            mut.append(total)
    if not wt or not mut:
        raise ValueError("Rosetta output missing WT or MUT rows")
    return (float(np.mean(mut)) - float(np.mean(wt))) / reu_per_kcal


# --------------------------------------------------------------------------- #
# Backends (adapter): real binaries on SCC vs. a local mock                      #
# --------------------------------------------------------------------------- #
class SubprocessBackend:
    """Runs the real FoldX / Rosetta binaries. Intended for the SCC only."""

    def run(self, method: str, record: dict, cfg: dict) -> list[float]:  # pragma: no cover
        raise NotImplementedError(
            "SubprocessBackend requires FoldX/Rosetta modules on the SCC. "
            "Set PRESCREEN_MOCK=1 to exercise the pipeline locally."
        )


class MockBackend:
    """Deterministic synthetic ΔΔG so the pipeline runs without the binaries.

    Values are reproducible per (method, variant) and are NOT science -- they exist
    only to exercise parsing/aggregation/IO. Never use mock output for a result.
    """

    def run(self, method: str, record: dict, cfg: dict) -> list[float]:
        n = cfg["prescreen"]["runs_per_variant"] if method == "foldx" else 1
        # stable_seed, not hash: hash is randomized per process, so "deterministic"
        # would have meant "different in every job" (see src/seeds.py).
        rng = np.random.default_rng(stable_seed(method, record["variant"]))
        center = 0.5 * (int(record["mature_pos"]) % 7) - 1.0  # arbitrary but deterministic
        return [float(center + rng.normal(0, 0.3)) for _ in range(n)]


def get_backend(cfg: dict):
    """MockBackend when ``PRESCREEN_MOCK`` is set, else the real SubprocessBackend."""
    return MockBackend() if os.environ.get("PRESCREEN_MOCK") else SubprocessBackend()


# --------------------------------------------------------------------------- #
# Per-variant prescreen                                                          #
# --------------------------------------------------------------------------- #
def _aggregate(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    return {
        "ddg": float(arr.mean()),
        "err": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
        "n": int(arr.size),
    }


def prescreen_variant(cfg: dict, variant: str, out_path: str | Path, backend=None) -> dict:
    """Run all configured methods for one variant and write ``results/prescreen/<v>.json``."""
    backend = backend or get_backend(cfg)
    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    methods = {}
    for method in cfg["prescreen"]["methods"]:
        methods[method] = _aggregate(backend.run(method, record, cfg))
    consensus = float(np.mean([m["ddg"] for m in methods.values()]))
    result = {"variant": variant, "mature_pos": int(record["mature_pos"]),
              "bucket": record["bucket"], "methods": methods, "consensus_ddg": consensus}

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))
    return result


# --------------------------------------------------------------------------- #
# Controls cross-check (advisory)                                               #
# --------------------------------------------------------------------------- #
def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    xr = np.argsort(np.argsort(x))
    yr = np.argsort(np.argsort(y))
    return _pearson(xr.astype(float), yr.astype(float))


def correlate_vs_experimental(results: list[dict], variants_csv: str | Path,
                              method: str = "consensus") -> dict:
    """Correlate prescreen ΔΔG against experimental ΔΔG on the controls.

    Advisory only -- the go/no-go gate is FEP's (README §2.2). Returns pearson,
    spearman, and n over controls that have a numeric ``exp_ddg``.
    """
    import csv
    exp = {}
    with open(variants_csv) as f:
        for row in csv.DictReader(f):
            if row["bucket"] == "positive_control" and row["exp_ddg"]:
                exp[row["variant"]] = float(row["exp_ddg"])
    pred, obs = [], []
    for r in results:
        if r["variant"] in exp:
            pred.append(r["consensus_ddg"] if method == "consensus" else r["methods"][method]["ddg"])
            obs.append(exp[r["variant"]])
    if len(pred) < 3:
        return {"pearson": None, "spearman": None, "n": len(pred)}
    x, y = np.asarray(pred), np.asarray(obs)
    return {"pearson": _pearson(x, y), "spearman": _spearman(x, y), "n": len(pred)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Empirical prescreen for one variant (Stage 2).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    result = prescreen_variant(cfg, args.variant, args.out)
    print(f"Prescreen {args.variant}: consensus ddG = {result['consensus_ddg']:.2f} kcal/mol -> {args.out}")


if __name__ == "__main__":
    main()
