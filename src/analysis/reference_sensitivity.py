"""CPU-only forensic sensitivity analysis for the failed SOD1 validation gate.

This module never writes a gate report and never changes ``data/variants.csv``. It freezes
the gate-attempt-1 predictions (including A4V's surviving run-time record), applies the
reference scenarios in ``data/reference_sensitivity_scenarios.csv``, and reports how the
descriptive metrics move. The output is not a second gate evaluation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from src.analysis.validate import accuracy_metrics, evaluate_gate
from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]

SCENARIOS = (
    ("gate_of_record", "Gate of record"),
    ("tavg_constant_dcp_half", "Dimer halved only (49.4 C)"),
    ("temp25_constant_dcp_dimer", "25 C only, constant dCp (dimer)"),
    ("temp25_constant_dcp_half", "25 C + halved, constant dCp"),
    ("tavg_tempdep_dcp_dimer", "49.4 C, temperature-dependent dCp (dimer)"),
    ("tavg_tempdep_dcp_half", "49.4 C + halved, temperature-dependent dCp"),
    ("temp25_tempdep_dcp_dimer", "25 C only, temperature-dependent dCp (dimer)"),
    ("temp25_tempdep_dcp_half", "25 C + halved, temperature-dependent dCp"),
)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Return Spearman correlation for the no-tie arrays used in this analysis."""
    xr = np.argsort(np.argsort(x)).astype(float)
    yr = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(xr, yr)[0, 1])


def sensitivity_metrics(pred: np.ndarray, obs: np.ndarray) -> dict[str, float]:
    """Return descriptive metrics without applying or changing gate thresholds."""
    metrics = accuracy_metrics(pred, obs)
    metrics.update({
        "pearson": float(np.corrcoef(pred, obs)[0, 1]),
        "spearman": _spearman(pred, obs),
        "sse": float(np.sum((pred - obs) ** 2)),
    })
    return metrics


def load_scenarios(path: str | Path) -> dict[str, dict[str, float]]:
    """Load the audited reference scenarios keyed by variant."""
    with open(path) as stream:
        rows = list(csv.DictReader(stream))
    required = {name for name, _ in SCENARIOS}
    if not rows or not required.issubset(rows[0]):
        missing = sorted(required - (set(rows[0]) if rows else set()))
        raise ValueError(f"reference scenario table is missing columns: {missing}")
    return {
        row["variant"]: {name: float(row[name]) for name in required}
        for row in rows
    }


def load_gate_attempt_1(archive_root: str | Path, gate_subset: list[str]) -> dict[str, dict]:
    """Load frozen attempt-1 records, using A4V's attested pre-resubmission JSON."""
    archive_root = Path(archive_root)
    records = {}
    for variant in gate_subset:
        path = (archive_root / "archive" / "A4V_3ns_ddg.json" if variant == "A4V"
                else archive_root / "fep" / variant / "ddg.json")
        if not path.exists():
            raise FileNotFoundError(f"missing frozen result for {variant}: {path}")
        record = json.loads(path.read_text())
        if record.get("variant") != variant:
            raise ValueError(f"{path} identifies {record.get('variant')!r}, expected {variant!r}")
        records[variant] = record
    return records


def analyze(archive_root: str | Path, scenario_path: str | Path,
            config_path: str | Path) -> tuple[list[str], dict[str, dict[str, float]], dict]:
    """Freeze the original usable subset and calculate all reference sensitivities."""
    cfg = load_config(config_path)
    subset = list(cfg["validation"]["gate_subset"])
    records = load_gate_attempt_1(archive_root, subset)
    refs = load_scenarios(scenario_path)

    missing = sorted(set(subset) - set(refs))
    if missing:
        raise ValueError(f"reference scenario table is missing gate variants: {missing}")

    original_exp = {variant: refs[variant]["gate_of_record"] for variant in subset}
    original_gate = evaluate_gate(records, original_exp, cfg)
    used = list(original_gate["used"])
    if used != [variant for variant in subset if variant != "F64A"]:
        raise ValueError(f"original exclusion set changed unexpectedly: used={used}")
    if original_gate["passed"] or original_gate["n"] != 7:
        raise ValueError("frozen inputs no longer reproduce the failed seven-point gate")

    pred = np.asarray([float(records[variant]["ddg"]) for variant in used])
    results = {}
    for name, _label in SCENARIOS:
        obs = np.asarray([refs[variant][name] for variant in used])
        results[name] = sensitivity_metrics(pred, obs)
    return used, results, original_gate


def plot_sensitivity(used: list[str], results: dict[str, dict[str, float]],
                     records: dict[str, dict], refs: dict[str, dict[str, float]],
                     output: str | Path) -> None:
    """Write the four-panel diagnostic comparison used by the methods deliverable."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    panels = (
        ("gate_of_record", "Gate of record"),
        ("tavg_constant_dcp_half", "Dimer normalization only"),
        ("temp25_constant_dcp_dimer", "Temperature only (constant dCp)"),
        ("temp25_constant_dcp_half", "Temperature + normalization"),
    )
    pred = np.asarray([float(records[variant]["ddg"]) for variant in used])
    colors = ["#b44747" if variant.startswith("G93") else "#31688e" for variant in used]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.2), sharex=True, sharey=True)
    for ax, (name, title) in zip(axes.flat, panels):
        obs = np.asarray([refs[variant][name] for variant in used])
        ax.scatter(obs, pred, c=colors, s=48, edgecolor="white", linewidth=0.7, zorder=3)
        for variant, x, y in zip(used, obs, pred):
            ax.annotate(variant, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
        ax.plot([-0.3, 7.3], [-0.3, 7.3], color="#555555", linestyle="--", linewidth=1)
        metric = results[name]
        ax.set_title(f"{title}\nr={metric['pearson']:.3f}; RMSE={metric['rmse']:.3f}")
        ax.grid(alpha=0.2)
        ax.set_xlim(-0.3, 7.3)
        ax.set_ylim(-0.3, 7.3)
    axes.flat[0].legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b44747",
               markeredgecolor="white", label="position 93"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#31688e",
               markeredgecolor="white", label="other sites"),
    ], loc="upper left", frameon=False, fontsize=8)
    fig.supxlabel("Experimental/reference delta-delta G (kcal/mol)", y=0.045)
    fig.supylabel("Frozen FEP prediction (kcal/mol)")
    fig.suptitle("Forensic reference sensitivity - not a re-evaluated validation gate", y=0.995)
    fig.text(0.5, 0.012,
             "G93S/G93V remain dimer-derived in every panel; halving is a normalization, "
             "not an isolated-monomer measurement.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0.02, 0.075, 1, 0.97))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", required=True,
                        help="root containing fep/ plus archive/A4V_3ns_ddg.json")
    parser.add_argument("--scenarios", default=ROOT / "data/reference_sensitivity_scenarios.csv",
                        type=Path)
    parser.add_argument("--config", default=ROOT / "config/pipeline.yaml", type=Path)
    parser.add_argument("--figure", type=Path, help="optional diagnostic figure output")
    args = parser.parse_args()

    used, results, original_gate = analyze(args.archive_root, args.scenarios, args.config)
    print(f"Frozen usable set: {', '.join(used)}")
    print("F64A exclusion preserved; original gate verdict preserved: "
          f"r={original_gate['pearson']:.6f}, RMSE={original_gate['rmse']:.6f}, FAILED")
    print("\n| scenario | Pearson r | Spearman rho | RMSE | MUE |")
    print("|---|---:|---:|---:|---:|")
    for name, label in SCENARIOS:
        metric = results[name]
        print(f"| {label} | {metric['pearson']:.3f} | {metric['spearman']:.3f} | "
              f"{metric['rmse']:.3f} | {metric['mue']:.3f} |")

    if args.figure:
        cfg = load_config(args.config)
        records = load_gate_attempt_1(args.archive_root, list(cfg["validation"]["gate_subset"]))
        refs = load_scenarios(args.scenarios)
        plot_sensitivity(used, results, records, refs, args.figure)


if __name__ == "__main__":
    main()
