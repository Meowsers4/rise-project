"""Stage 5 -- the validation GATE, then classification (README §2.2, CLAUDE.md #2).

Go/no-go: FEP ΔΔG must reproduce the experimental ΔΔG of the control ``gate_subset``
above ``validation.min_pearson`` BEFORE any uncharacterized variant is trusted. The
threshold is read from config and never lowered here. If the gate fails, this module
writes the (still presentable) control-correlation report and refuses to emit a
classification map -- extending to novel variants on a failed gate is confident
nonsense.

FEP per-variant results are read from ``results/fep/<variant>/ddg.json`` with the
schema :data:`FEP_DDG_SCHEMA` (written by ``src/fep/analyze.py``):

    {"variant": str, "ddg": float, "ddg_err": float,
     "cycle_closure_kcal": float, "converged": bool, "provenance": str}

``provenance`` must equal the configured ``fep.framework``. Mock output, and any file
that carries no provenance at all (i.e. one written or edited by hand), is excluded from
the gate and reported in ``excluded`` -- a gate is only as trustworthy as the weakest
number in it, and a plausible-looking ΔΔG of unknown origin is the failure mode this
whole module exists to prevent.

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

FEP_DDG_SCHEMA = ("variant", "ddg", "ddg_err", "cycle_closure_kcal", "converged",
                  "provenance")


def is_production_result(rec: dict, cfg: dict) -> bool:
    """True if ``rec`` came from the configured production engine (``fep.framework``).

    A missing ``provenance`` is never production: that is what an unlabelled, mock, or
    hand-written ΔΔG looks like, and it must not reach the gate.
    """
    return rec.get("provenance") == cfg["fep"]["framework"]


# --------------------------------------------------------------------------- #
# Small stats (numpy-only)                                                      #
# --------------------------------------------------------------------------- #
def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    xr = np.argsort(np.argsort(x)).astype(float)
    yr = np.argsort(np.argsort(y)).astype(float)
    return _pearson(xr, yr)


def accuracy_metrics(pred: np.ndarray, obs: np.ndarray) -> dict:
    """Absolute-accuracy metrics to sit alongside the correlation.

    Correlation answers "do these lie on *a* line?", not "on the *right* line". A method
    reporting exactly half of every true ΔΔG scores r = 1.0 while every number is 50%
    wrong. RMSE/MUE are the "can I trust one prediction" numbers; slope and intercept say
    whether the errors are a systematic scaling/offset rather than scatter.
    """
    err = pred - obs
    slope, intercept = np.polyfit(obs, pred, 1) if len(obs) > 1 else (float("nan"),) * 2
    return {
        "rmse": float(np.sqrt(np.mean(err ** 2))),
        "mue": float(np.mean(np.abs(err))),
        "max_abs_err": float(np.max(np.abs(err))),
        "slope": float(slope),
        "intercept": float(intercept),
    }


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

    A gate variant is *usable* only if it came from the production engine
    (:func:`is_production_result`), converged, and has cycle closure within
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
        if not is_production_result(rec, cfg):
            excluded.append({
                "variant": v,
                "reason": f"provenance {rec.get('provenance', 'MISSING')!r} is not the "
                          f"production engine {cfg['fep']['framework']!r}",
            })
            continue
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

    pred_a, obs_a = np.asarray(pred), np.asarray(obs)
    pearson = _pearson(pred_a, obs_a)
    spearman = _spearman(pred_a, obs_a)
    acc = accuracy_metrics(pred_a, obs_a)

    max_rmse = vcfg.get("max_rmse_kcal")          # optional second criterion
    checks = {"pearson": pearson >= min_pearson}
    if max_rmse is not None:
        checks["rmse"] = acc["rmse"] <= max_rmse
    passed = all(checks.values())
    failed = [k for k, ok in checks.items() if not ok]
    reason = "correlation and accuracy criteria met" if passed else \
             f"failed on: {', '.join(failed)}"

    return {"passed": bool(passed), "reason": reason,
            "pearson": pearson, "spearman": spearman, **acc, "n": n,
            "min_pearson": min_pearson, "max_rmse_kcal": max_rmse,
            "used": used, "excluded": excluded}


def classify_uncharacterized(fep: dict[str, dict], panel: list[dict], cfg: dict) -> list[dict]:
    """Label uncharacterized variants destabilizing/stable by FEP ΔΔG vs threshold.

    Non-production results are skipped for the same reason the gate excludes them: a
    classification is a claim about a variant, and a mock ΔΔG cannot support one.
    """
    threshold = cfg["validation"]["destabilizing_ddg_kcal"]
    rows = []
    for row in panel:
        v = row["variant"]
        if row["bucket"] == "positive_control" or v not in fep:
            continue
        if not is_production_result(fep[v], cfg):
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
            f"rmse={gate.get('rmse')}, n={gate['n']}). Refusing to classify "
            f"uncharacterized variants. See {out_gate_report}."
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
    print(f"GATE PASSED (n={gate['n']}): pearson={gate['pearson']:.3f} "
          f"(min {gate['min_pearson']}), RMSE={gate['rmse']:.2f} kcal/mol "
          f"(max {gate['max_rmse_kcal']}), MUE={gate['mue']:.2f}, "
          f"slope={gate['slope']:.2f}. Wrote {args.out}.")


if __name__ == "__main__":
    main()
