"""Stage 5 -- the validation GATE, then classification (README §2.2, CLAUDE.md #2).

Go/no-go: FEP ΔΔG must reproduce the experimental ΔΔG of the control ``gate_subset``
above ``validation.min_pearson`` BEFORE any uncharacterized variant is trusted. The
threshold is read from config and never lowered here. If the gate fails, this module
writes the (still presentable) control-correlation report and refuses to emit a
classification map -- extending to novel variants on a failed gate is confident
nonsense.

FEP per-variant results are read from ``results/fep/<variant>/ddg.json`` with the
schema :data:`FEP_DDG_SCHEMA` (the contract for the future ``src/fep/analyze.py``):

    {"variant": str, "ddg": float, "ddg_err": float,
     "cycle_closure_kcal": float, "converged": bool}

ΔΔG sign is positive = destabilizing everywhere (``validation.ddg_sign``).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]

FEP_DDG_SCHEMA = ("variant", "ddg", "ddg_err", "cycle_closure_kcal", "converged")


# --------------------------------------------------------------------------- #
# Small stats (numpy-only)                                                      #
# --------------------------------------------------------------------------- #
def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    xr = np.argsort(np.argsort(x)).astype(float)
    yr = np.argsort(np.argsort(y)).astype(float)
    return _pearson(xr, yr)


# --------------------------------------------------------------------------- #
# IO                                                                            #
# --------------------------------------------------------------------------- #
def load_fep_results(fep_dir: str | Path, variants: list[str]) -> dict[str, dict]:
    """Load ``results/fep/<variant>/ddg.json`` for each variant that has one."""
    fep_dir = Path(fep_dir)
    out: dict[str, dict] = {}
    for v in variants:
        path = fep_dir / v / "ddg.json"
        if path.exists():
            out[v] = json.loads(path.read_text())
    return out


def load_panel(variants_csv: str | Path) -> list[dict]:
    with open(variants_csv) as f:
        return list(csv.DictReader(f))


# --------------------------------------------------------------------------- #
# The gate                                                                      #
# --------------------------------------------------------------------------- #
def evaluate_gate(fep: dict[str, dict], exp: dict[str, float], cfg: dict) -> dict:
    """Correlate FEP vs experimental ΔΔG on the gate_subset; return go/no-go.

    A gate variant is *usable* only if it converged and its cycle closure is within
    ``fep.convergence.max_cycle_closure_kcal``. The gate never lowers min_pearson.
    """
    vcfg = cfg["validation"]
    subset = vcfg.get("gate_subset") or list(exp.keys())
    max_closure = cfg["fep"]["convergence"]["max_cycle_closure_kcal"]
    min_pearson = vcfg["min_pearson"]
    min_points = vcfg.get("min_gate_points", 3)

    used, excluded = [], []
    pred, obs = [], []
    for v in subset:
        if v not in fep or v not in exp:
            excluded.append({"variant": v, "reason": "missing FEP or experimental value"})
            continue
        rec = fep[v]
        if not rec.get("converged", False) or rec.get("cycle_closure_kcal", 1e9) > max_closure:
            excluded.append({"variant": v, "reason": "not converged / cycle closure too large"})
            continue
        used.append(v)
        pred.append(rec["ddg"])
        obs.append(exp[v])

    n = len(used)
    if n < min_points:
        return {"passed": False, "reason": f"only {n} usable gate points (< {min_points})",
                "pearson": None, "spearman": None, "n": n,
                "min_pearson": min_pearson, "used": used, "excluded": excluded}

    pearson = _pearson(np.asarray(pred), np.asarray(obs))
    spearman = _spearman(np.asarray(pred), np.asarray(obs))
    passed = pearson >= min_pearson
    return {"passed": bool(passed),
            "reason": "pearson >= min_pearson" if passed else "pearson below min_pearson",
            "pearson": pearson, "spearman": spearman, "n": n,
            "min_pearson": min_pearson, "used": used, "excluded": excluded}


def classify_uncharacterized(fep: dict[str, dict], panel: list[dict], cfg: dict) -> list[dict]:
    """Label uncharacterized variants destabilizing/stable by FEP ΔΔG vs threshold."""
    threshold = cfg["validation"]["destabilizing_ddg_kcal"]
    rows = []
    for row in panel:
        v = row["variant"]
        if row["bucket"] == "positive_control" or v not in fep:
            continue
        ddg = fep[v]["ddg"]
        rows.append({
            "variant": v,
            "bucket": row["bucket"],
            "oligomer": row["oligomer"],
            "fep_ddg": ddg,
            "fep_ddg_err": fep[v].get("ddg_err"),
            "prediction": "destabilizing" if ddg > threshold else "stable",
        })
    return rows


def run_validation(cfg: dict, out_ddg_map: str | Path, out_gate_report: str | Path,
                   fep_dir: str | Path | None = None,
                   variants_csv: str | Path | None = None) -> dict:
    """Run the gate; write the report always and the ddG map only if the gate passes.

    Returns the gate report dict. Raises SystemExit on no-go so the Snakemake rule
    fails and the pipeline stops before trusting uncharacterized variants.
    """
    variants_csv = variants_csv or ROOT / cfg["panel"]["csv"]
    fep_dir = fep_dir or ROOT / "results" / "fep"
    panel = load_panel(variants_csv)
    variants = [r["variant"] for r in panel]
    exp = {r["variant"]: float(r["exp_ddg"])
           for r in panel if r["bucket"] == "positive_control" and r["exp_ddg"]}
    fep = load_fep_results(fep_dir, variants)

    gate = evaluate_gate(fep, exp, cfg)

    out_gate_report = Path(out_gate_report)
    out_gate_report.parent.mkdir(parents=True, exist_ok=True)
    out_gate_report.write_text(json.dumps(gate, indent=2))

    if not gate["passed"]:
        raise SystemExit(
            f"VALIDATION GATE FAILED ({gate['reason']}; pearson={gate['pearson']}, "
            f"n={gate['n']}). Refusing to classify uncharacterized variants. "
            f"See {out_gate_report}."
        )

    classifications = classify_uncharacterized(fep, panel, cfg)
    out_ddg_map = Path(out_ddg_map)
    out_ddg_map.parent.mkdir(parents=True, exist_ok=True)
    with open(out_ddg_map, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["variant", "bucket", "oligomer", "fep_ddg", "fep_ddg_err", "prediction"]
        )
        writer.writeheader()
        writer.writerows(classifications)
    return gate


def main() -> None:
    parser = argparse.ArgumentParser(description="Validation gate + classification (Stage 5).")
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True, help="ddg_map.csv path")
    args = parser.parse_args()

    cfg = load_config(args.config)
    gate = run_validation(cfg, args.out, ROOT / cfg["validation"]["outputs"]["gate_report"])
    print(f"GATE PASSED: pearson={gate['pearson']:.3f} (n={gate['n']}, "
          f"min={gate['min_pearson']}). Wrote {args.out}.")


if __name__ == "__main__":
    main()
