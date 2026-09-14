"""Build the CPU-only forensic figure package from archived SOD1 results.

The builders read the frozen gate-attempt-1 records, independently re-derived
convergence JSON, the A4V 9 ns diagnostic, and the G93A SS archive.  They never
write a gate report or modify a simulation result.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from src.analysis.reference_sensitivity import (
    load_gate_attempt_1,
    load_scenarios,
)
from src.analysis.validate import evaluate_gate
from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]


def _read_json(path: str | Path) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def load_convergence_records(path: str | Path) -> dict[str, dict]:
    """Load one convergence JSON per variant from ``path``."""
    path = Path(path)
    records: dict[str, dict] = {}
    for json_path in sorted(path.glob("*.json")):
        record = _read_json(json_path)
        variant = record.get("variant")
        if not variant:
            raise ValueError(f"missing variant in {json_path}")
        if variant in records:
            raise ValueError(f"duplicate convergence record for {variant}")
        records[variant] = record
    if not records:
        raise ValueError(f"no convergence JSON files found in {path}")
    return records


def folded_disagreement_rows(records: dict[str, dict]) -> list[dict]:
    """Return folded hysteresis and independent-box disagreement per replicate.

    Disagreement is the definition frozen in ``stop_rule_reanalysis.md``:
    ``abs(DG_i - mean(DG_j for j != i))`` within a variant and leg.
    """
    rows: list[dict] = []
    for variant, record in sorted(records.items()):
        legs = sorted(
            (leg for leg in record.get("legs", []) if leg.get("leg") == "folded"),
            key=lambda leg: int(leg["replicate"]),
        )
        if len(legs) != 3:
            raise ValueError(f"{variant} has {len(legs)} folded replicates, expected 3")
        values = np.asarray([float(leg["dg_kcal"]) for leg in legs])
        for index, leg in enumerate(legs):
            siblings = np.delete(values, index)
            rows.append({
                "variant": variant,
                "replicate": int(leg["replicate"]),
                "dg_kcal": float(leg["dg_kcal"]),
                "hysteresis_kcal": float(leg["hysteresis_kcal"]),
                "disagreement_kcal": float(abs(values[index] - np.mean(siblings))),
            })
    return rows


def disagreement_correlations(rows: list[dict]) -> tuple[float, float]:
    """Return Pearson and no-tie Spearman correlations for diagnostic rows."""
    hysteresis = np.asarray([row["hysteresis_kcal"] for row in rows])
    disagreement = np.asarray([row["disagreement_kcal"] for row in rows])
    pearson = float(np.corrcoef(hysteresis, disagreement)[0, 1])
    hx = np.argsort(np.argsort(hysteresis)).astype(float)
    dy = np.argsort(np.argsort(disagreement)).astype(float)
    spearman = float(np.corrcoef(hx, dy)[0, 1])
    return pearson, spearman


def _leg_rows(record: dict, leg_name: str) -> list[dict]:
    rows = sorted(
        (leg for leg in record.get("legs", []) if leg.get("leg") == leg_name),
        key=lambda leg: int(leg["replicate"]),
    )
    if len(rows) != 3:
        raise ValueError(
            f"{record.get('variant', '<unknown>')} has {len(rows)} {leg_name} replicates, "
            "expected 3"
        )
    return rows


def inspect_cysteine_topology(path: str | Path) -> dict:
    """Report every cysteine SG--SG bond and cysteine retaining an HG atom."""
    section = ""
    atoms: dict[int, tuple[int, str, str]] = {}
    bonds: set[frozenset[int]] = set()
    for raw_line in Path(path).read_text().splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            continue
        fields = line.split()
        if section == "atoms" and len(fields) >= 5:
            try:
                atom_index = int(fields[0])
                residue_number = int(fields[2])
            except ValueError:
                continue
            atoms[atom_index] = (residue_number, fields[3], fields[4])
        elif section == "bonds" and len(fields) >= 2:
            try:
                bonds.add(frozenset((int(fields[0]), int(fields[1]))))
            except ValueError:
                continue

    sg_by_residue: dict[int, int] = {}
    hg_residues: set[int] = set()
    for atom_index, (residue_number, residue_name, atom_name) in atoms.items():
        if residue_name != "CYS":
            continue
        if atom_name == "SG":
            sg_by_residue[residue_number] = atom_index
        elif atom_name == "HG":
            hg_residues.add(residue_number)

    residue_by_sg = {atom_index: residue for residue, atom_index in sg_by_residue.items()}
    sg_bonds = sorted(
        tuple(sorted((residue_by_sg[first], residue_by_sg[second])))
        for bond in bonds if len(bond) == 2
        for first, second in [tuple(bond)]
        if first in residue_by_sg and second in residue_by_sg
    )
    return {
        "sg_bonds": sg_bonds,
        "hg_residues": sorted(hg_residues),
        "sg_atom_indices": sg_by_residue,
    }


def parse_cys_topology(path: str | Path, residues: tuple[int, int] = (57, 146)) -> dict:
    """Inspect one GROMACS topology for the Cys57--Cys146 bond and thiol H atoms."""
    topology = inspect_cysteine_topology(path)
    sg_by_residue = topology["sg_atom_indices"]
    if not set(residues).issubset(sg_by_residue):
        raise ValueError(f"could not find both target SG atoms in {path}")
    target_bond = tuple(sorted(residues))
    return {
        "bonded": target_bond in topology["sg_bonds"],
        "hg_present": {
            residue: residue in topology["hg_residues"] for residue in residues
        },
        "sg_atom_indices": {residue: sg_by_residue[residue] for residue in residues},
    }


def inspect_topology_set(root: str | Path) -> list[dict]:
    """Inspect all three folded replicate topologies beneath an archived FEP tree."""
    root = Path(root)
    paths = sorted(root.glob("folded/system_r*/hybrid.top"))
    if len(paths) != 3:
        raise ValueError(f"found {len(paths)} folded hybrid.top files under {root}, expected 3")
    return [{"path": path, **parse_cys_topology(path)} for path in paths]


def g93a_redox_summary(baseline_record: dict, ss_record: dict) -> dict:
    """Return the primary folded-leg comparison for the G93A redox diagnostic."""
    baseline = np.asarray([float(row["dg_kcal"]) for row in _leg_rows(baseline_record, "folded")])
    ss = np.asarray([float(row["dg_kcal"]) for row in _leg_rows(ss_record, "folded")])
    return {
        "baseline": baseline,
        "ss": ss,
        "baseline_mean": float(np.mean(baseline)),
        "ss_mean": float(np.mean(ss)),
        "shift": float(np.mean(ss) - np.mean(baseline)),
    }


def _style_axes(ax) -> None:
    ax.grid(color="#d8d8d8", linewidth=0.7, alpha=0.65, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_gate_of_record(records: dict[str, dict], refs: dict[str, dict[str, float]],
                        gate: dict, output: str | Path) -> None:
    """Plot the immutable gate of record, including the pre-registered F64A exclusion."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    used = set(gate["used"])
    variants = list(records)
    fig, ax = plt.subplots(figsize=(7.5, 6.7))
    for variant in variants:
        x = float(refs[variant]["gate_of_record"])
        y = float(records[variant]["ddg"])
        err = float(records[variant].get("ddg_err", 0.0))
        if variant in used:
            color = "#b44747" if variant.startswith("G93") else "#31688e"
            ax.errorbar(x, y, yerr=err, fmt="o", ms=7, color=color,
                        ecolor=color, capsize=2.5, linewidth=1.1, zorder=3)
        else:
            ax.errorbar(x, y, yerr=err, fmt="X", ms=9, color="#555555",
                        ecolor="#777777", capsize=2.5, linewidth=1.1, zorder=3)
        offset = (7, -13) if variant == "F64A" else (5, 5)
        ax.annotate(variant, (x, y), xytext=offset, textcoords="offset points", fontsize=9)

    low = min(min(float(refs[v]["gate_of_record"]) for v in variants),
              min(float(records[v]["ddg"]) for v in variants)) - 0.25
    high = max(max(float(refs[v]["gate_of_record"]) for v in variants),
               max(float(records[v]["ddg"]) for v in variants)) + 0.35
    ax.plot([low, high], [low, high], color="#555555", linestyle="--", linewidth=1)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Experimental reference ΔΔG (kcal/mol)")
    ax.set_ylabel("FEP prediction ΔΔG (kcal/mol)")
    ax.set_title("Pre-registered validation gate of record — FAILED")
    ax.text(0.98, 0.04,
            f"usable n = {gate['n']}\nPearson r = {gate['pearson']:.3f}\n"
            f"RMSE = {gate['rmse']:.3f} kcal/mol",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "pad": 5})
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b44747",
               markeredgecolor="#b44747", label="usable, position 93"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#31688e",
               markeredgecolor="#31688e", label="usable, other sites"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor="#555555",
               markeredgecolor="#555555", label="F64A excluded by closure rule"),
    ], loc="center right", frameon=False, fontsize=8.5)
    _style_axes(ax)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_convergence_diagnostic(rows: list[dict], output: str | Path) -> None:
    """Plot folded hysteresis against independently solvated-box disagreement."""
    import matplotlib.pyplot as plt

    pearson, spearman = disagreement_correlations(rows)
    fig, ax = plt.subplots(figsize=(7.6, 6.3))
    for row in rows:
        is_f64a = row["variant"] == "F64A"
        ax.scatter(row["hysteresis_kcal"], row["disagreement_kcal"],
                   s=66 if is_f64a else 42, marker="^" if is_f64a else "o",
                   color="#b44747" if is_f64a else "#31688e",
                   alpha=0.95 if is_f64a else 0.75, edgecolor="white", linewidth=0.6,
                   zorder=3)
        if is_f64a and row["replicate"] == 1:
            ax.annotate("F64A r1\nrank 1/48 disagreement",
                        (row["hysteresis_kcal"], row["disagreement_kcal"]),
                        xytext=(28, -8), textcoords="offset points", fontsize=9,
                        arrowprops={"arrowstyle": "-", "color": "#777777"})
    ax.set_xlabel("Within-ladder hysteresis (kcal/mol)")
    ax.set_ylabel("Independent-box disagreement (kcal/mol)")
    ax.set_title("Folded-leg diagnostics do not separate replicate disagreement")
    ax.text(0.98, 0.95,
            f"n = {len(rows)}\nPearson r = {pearson:+.3f}\nSpearman ρ = {spearman:+.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "pad": 5})
    ax.scatter([], [], marker="o", color="#31688e", label="other records")
    ax.scatter([], [], marker="^", color="#b44747", label="F64A records")
    ax.legend(frameon=False, loc="center right", fontsize=8.5)
    _style_axes(ax)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_a4v_sampling(attempt1: dict, attempt2: dict, experimental: float,
                      output: str | Path) -> None:
    """Plot the A4V 3 ns versus 9 ns precision and accuracy diagnostics."""
    import matplotlib.pyplot as plt

    labels = ("3 ns", "9 ns")
    records = (attempt1, attempt2)
    colors = ("#8c6bb1", "#238b8e")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.9))

    folded_values = []
    for index, (record, color) in enumerate(zip(records, colors)):
        legs = _leg_rows(record, "folded")
        values = np.asarray([float(row["dg_kcal"]) for row in legs])
        errors = np.asarray([float(row["dg_err_kcal"]) for row in legs])
        folded_values.append(values)
        jitter = np.asarray([-0.055, 0.0, 0.055])
        axes[0].errorbar(index + jitter, values, yerr=errors, fmt="o", color=color,
                         ecolor=color, capsize=2.5, linewidth=1, ms=6, zorder=3)
        axes[0].hlines(np.mean(values), index - 0.18, index + 0.18,
                       color=color, linewidth=2.3, zorder=2)
        axes[0].text(index, min(values) - 0.18,
                     f"spread {np.ptp(values):.2f}", ha="center", va="top", fontsize=9)
    for replicate in range(3):
        axes[0].plot([0 - 0.055 + 0.055 * replicate, 1 - 0.055 + 0.055 * replicate],
                     [folded_values[0][replicate], folded_values[1][replicate]],
                     color="#bdbdbd", linewidth=0.8, zorder=1)
    folded_all = np.concatenate(folded_values)
    axes[0].set_ylim(float(np.min(folded_all)) - 0.28,
                     float(np.max(folded_all)) + 0.17)
    axes[0].set_xticks((0, 1), labels)
    axes[0].set_ylabel("Folded ΔG (kcal/mol)")
    axes[0].set_title("Independent-box precision")
    _style_axes(axes[0])

    ddg_values = [np.asarray(record["per_replicate_ddg"], dtype=float) for record in records]
    for replicate in range(3):
        axes[1].plot((0, 1), (ddg_values[0][replicate], ddg_values[1][replicate]),
                     color="#bdbdbd", linewidth=0.8, zorder=1)
    for index, (record, values, color) in enumerate(zip(records, ddg_values, colors)):
        axes[1].scatter(np.full(3, index) + np.asarray([-0.055, 0, 0.055]), values,
                        color=color, s=42, edgecolor="white", linewidth=0.6, zorder=3)
        axes[1].errorbar(index, float(record["ddg"]), yerr=float(record["ddg_err"]),
                         fmt="D", color=color, ms=6, capsize=3, linewidth=1.2, zorder=4)
        error = float(record["ddg"]) - experimental
        axes[1].text(index, max(values) + 0.14, f"error {error:+.2f}",
                     ha="center", va="bottom", fontsize=9)
    axes[1].axhline(experimental, color="#555555", linestyle="--", linewidth=1.2,
                    label=f"experiment {experimental:.2f}")
    axes[1].set_ylim(experimental - 0.12,
                     float(max(np.max(values) for values in ddg_values)) + 0.36)
    axes[1].set_xticks((0, 1), labels)
    axes[1].set_ylabel("A4V ΔΔG (kcal/mol)")
    axes[1].set_title("Accuracy remains outside the reference")
    axes[1].legend(frameon=False, loc="lower right", fontsize=8.5)
    _style_axes(axes[1])

    fig.suptitle("A4V: tripled folded sampling improves precision, not accuracy", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_g93a_redox(summary: dict, baseline_topologies: list[dict],
                    ss_topologies: list[dict], output: str | Path) -> None:
    """Plot the negative G93A folded-state disulfide diagnostic."""
    import matplotlib.pyplot as plt

    baseline = summary["baseline"]
    ss = summary["ss"]
    if not all(not row["bonded"] and all(row["hg_present"].values())
               for row in baseline_topologies):
        raise ValueError("2SH topology set does not show two reduced target cysteines")
    if not all(row["bonded"] and not any(row["hg_present"].values())
               for row in ss_topologies):
        raise ValueError("SS topology set does not show the C57--C146 disulfide")

    fig, ax = plt.subplots(figsize=(7.3, 5.9))
    colors = ("#31688e", "#b44747")
    for index, (values, color) in enumerate(((baseline, colors[0]), (ss, colors[1]))):
        x = np.full(3, index) + np.asarray([-0.06, 0.0, 0.06])
        ax.scatter(x, values, s=55, color=color, edgecolor="white", linewidth=0.7, zorder=3)
        ax.errorbar(index, np.mean(values), yerr=np.std(values, ddof=1), fmt="D",
                    color=color, capsize=4, linewidth=1.3, ms=6, zorder=4)
    ax.set_xticks((0, 1), ("apo-2SH baseline", "apo-SS diagnostic"))
    ax.set_ylabel("G93A folded ΔG (kcal/mol)")
    ax.set_title("Restoring C57–C146 does not shift G93A's folded mutation cost")
    ax.text(0.5, 0.94,
            f"SS − 2SH mean = {summary['shift']:+.4f} kcal/mol",
            transform=ax.transAxes, ha="center", va="top", fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "pad": 5})
    ax.text(0.5, 0.05,
            "Topology gate (3/3 boxes each): 2SH has no SG–SG bond and both HG;\n"
            "SS has the C57–C146 SG–SG bond and neither HG.",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5)
    _style_axes(ax)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def build_figure_package(attempt1_archive: str | Path,
                         reanalysis_convergence: str | Path,
                         ss_archive: str | Path,
                         scenarios_path: str | Path,
                         config_path: str | Path,
                         output_dir: str | Path) -> dict:
    """Build Figures 1--3 and 5 and return their source-derived summaries."""
    attempt1_archive = Path(attempt1_archive)
    ss_archive = Path(ss_archive)
    output_dir = Path(output_dir)
    cfg = load_config(config_path)
    variants = list(cfg["validation"]["gate_subset"])
    records = load_gate_attempt_1(attempt1_archive, variants)
    refs = load_scenarios(scenarios_path)
    experimental = {variant: refs[variant]["gate_of_record"] for variant in variants}
    gate = evaluate_gate(records, experimental, cfg)
    if gate["passed"] or gate["n"] != 7:
        raise ValueError("attempt-1 archive no longer reproduces the failed seven-point gate")

    convergence = load_convergence_records(reanalysis_convergence)
    disagreement = folded_disagreement_rows(convergence)
    if len(disagreement) != 24:
        raise ValueError(f"expected 24 folded diagnostic rows, found {len(disagreement)}")
    pearson, spearman = disagreement_correlations(disagreement)

    a4v_3ns = _read_json(attempt1_archive / "archive" / "A4V_3ns_convergence.json")
    a4v_9ns = _read_json(attempt1_archive / "convergence" / "A4V.json")
    g93a_2sh = convergence["G93A"]
    g93a_ss = _read_json(ss_archive / "convergence_SS.json")
    redox = g93a_redox_summary(g93a_2sh, g93a_ss)
    baseline_topologies = inspect_topology_set(attempt1_archive / "fep" / "G93A")
    ss_topologies = inspect_topology_set(ss_archive)

    plot_gate_of_record(records, refs, gate, output_dir / "gate_of_record.png")
    plot_convergence_diagnostic(disagreement, output_dir / "convergence_diagnostic.png")
    plot_a4v_sampling(a4v_3ns, a4v_9ns, experimental["A4V"],
                      output_dir / "a4v_sampling_sensitivity.png")
    plot_g93a_redox(redox, baseline_topologies, ss_topologies,
                    output_dir / "g93a_disulfide_diagnostic.png")
    return {
        "gate": gate,
        "diagnostic_pearson": pearson,
        "diagnostic_spearman": spearman,
        "diagnostic_n": len(disagreement),
        "g93a_folded_shift_kcal": redox["shift"],
        "outputs": [
            output_dir / "gate_of_record.png",
            output_dir / "convergence_diagnostic.png",
            output_dir / "a4v_sampling_sensitivity.png",
            output_dir / "g93a_disulfide_diagnostic.png",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt1-archive", required=True, type=Path)
    parser.add_argument("--reanalysis-convergence", required=True, type=Path)
    parser.add_argument("--ss-archive", required=True, type=Path)
    parser.add_argument("--scenarios", default=ROOT / "data/reference_sensitivity_scenarios.csv",
                        type=Path)
    parser.add_argument("--config", default=ROOT / "config/pipeline.yaml", type=Path)
    parser.add_argument("--output-dir", default=ROOT / "docs/figures", type=Path)
    args = parser.parse_args()
    summary = build_figure_package(
        args.attempt1_archive,
        args.reanalysis_convergence,
        args.ss_archive,
        args.scenarios,
        args.config,
        args.output_dir,
    )
    print("Built CPU-only forensic figure package:")
    for path in summary["outputs"]:
        print(f"  {path}")
    print(
        f"Checks: gate r={summary['gate']['pearson']:.6f}, "
        f"RMSE={summary['gate']['rmse']:.6f}; "
        f"folded diagnostic r={summary['diagnostic_pearson']:+.6f}, "
        f"rho={summary['diagnostic_spearman']:+.6f}; "
        f"G93A SS-2SH={summary['g93a_folded_shift_kcal']:+.6f} kcal/mol"
    )


if __name__ == "__main__":
    main()
