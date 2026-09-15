#!/usr/bin/env python3
"""Independent reproduction of reference-sensitivity Table 3 and Figure 4.

This script intentionally imports none of ``src.analysis``.  It reads the frozen
run records and audited CSV directly, applies the documented per-variant closure
cap, and implements Pearson correlation, Spearman correlation, RMSE, and MUE with
the Python standard library.  Matplotlib is used only when ``--figure`` is given.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


GATE_VARIANTS = ("I18V", "I113T", "A4V", "G93A", "G93S", "F64A", "I149A", "G93V")
CLOSURE_CAP_KCAL = 1.0
SCENARIOS = (
    ("gate_of_record", "gate of record"),
    ("tavg_constant_dcp_half", "dimer normalization only, constant dCp"),
    ("temp25_constant_dcp_dimer", "temperature only, 25 C constant dCp, still dimer"),
    ("temp25_constant_dcp_half", "temperature + normalization, constant dCp"),
    ("tavg_tempdep_dcp_dimer", "49.4 C, temperature-dependent dCp, still dimer"),
    ("tavg_tempdep_dcp_half", "normalization only, temperature-dependent dCp"),
    ("temp25_tempdep_dcp_dimer", "temperature only, 25 C temperature-dependent dCp, still dimer"),
    ("temp25_tempdep_dcp_half", "temperature + normalization, temperature-dependent dCp"),
)


def _mean(values: list[float]) -> float:
    if not values:
        raise ValueError("mean requires at least one value")
    return sum(values) / len(values)


def pearson(x: list[float], y: list[float]) -> float:
    """Return the sample Pearson correlation from its defining sums."""
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("Pearson inputs must have the same length >= 2")
    x_mean, y_mean = _mean(x), _mean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    x_ss = sum((a - x_mean) ** 2 for a in x)
    y_ss = sum((b - y_mean) ** 2 for b in y)
    return numerator / math.sqrt(x_ss * y_ss)


def ranks(values: list[float]) -> list[float]:
    """Return average ranks, including a tie-safe path not needed by these data."""
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    result = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        stop = start + 1
        while stop < len(ordered) and ordered[stop][1] == ordered[start][1]:
            stop += 1
        average_rank = ((start + 1) + stop) / 2.0
        for position in range(start, stop):
            result[ordered[position][0]] = average_rank
        start = stop
    return result


def metrics(predicted: list[float], observed: list[float]) -> dict[str, float]:
    """Calculate the four quantities printed in the audited sensitivity table."""
    if len(predicted) != len(observed) or not predicted:
        raise ValueError("metric inputs must have equal non-zero length")
    errors = [pred - obs for pred, obs in zip(predicted, observed)]
    return {
        "pearson": pearson(predicted, observed),
        "spearman": pearson(ranks(predicted), ranks(observed)),
        "rmse": math.sqrt(_mean([error * error for error in errors])),
        "mue": _mean([abs(error) for error in errors]),
    }


def load_records(archive_root: Path) -> dict[str, dict]:
    """Load the eight frozen records without using project loader functions."""
    records: dict[str, dict] = {}
    for variant in GATE_VARIANTS:
        if variant == "A4V":
            path = archive_root / "archive" / "A4V_3ns_ddg.json"
        else:
            path = archive_root / "fep" / variant / "ddg.json"
        record = json.loads(path.read_text())
        if record.get("variant") != variant:
            raise ValueError(f"{path} identifies {record.get('variant')!r}, expected {variant}")
        records[variant] = record
    return records


def load_references(path: Path) -> dict[str, dict[str, float]]:
    """Load the audited scenario matrix with explicit numeric conversion."""
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required_columns = {"variant", *(name for name, _ in SCENARIOS)}
    actual_columns = set(rows[0]) if rows else set()
    if not rows or not required_columns.issubset(actual_columns):
        raise ValueError(f"missing scenario columns: {sorted(required_columns - actual_columns)}")
    references = {
        row["variant"]: {name: float(row[name]) for name, _ in SCENARIOS}
        for row in rows
    }
    if set(references) != set(GATE_VARIANTS):
        raise ValueError("scenario variants do not match the frozen eight-variant panel")
    return references


def reproduce(archive_root: Path, scenario_path: Path) -> dict:
    """Return the independently calculated usable set, inputs, and metrics."""
    records = load_records(archive_root)
    references = load_references(scenario_path)
    used = [
        variant for variant in GATE_VARIANTS
        if float(records[variant]["cycle_closure_kcal"]) <= CLOSURE_CAP_KCAL
    ]
    excluded = [variant for variant in GATE_VARIANTS if variant not in used]
    if excluded != ["F64A"]:
        raise ValueError(f"independent closure filter excluded {excluded}, expected only F64A")

    predicted = [float(records[variant]["ddg"]) for variant in used]
    results = {}
    for name, label in SCENARIOS:
        observed = [references[variant][name] for variant in used]
        results[name] = {"label": label, **metrics(predicted, observed)}
    return {
        "used": used,
        "excluded": excluded,
        "closure_cap_kcal": CLOSURE_CAP_KCAL,
        "predicted": dict(zip(used, predicted)),
        "metrics": results,
    }


def plot(reproduction: dict, scenario_path: Path, output: Path) -> None:
    """Draw an independent four-panel version of Figure 4."""
    import matplotlib.pyplot as plt

    references = load_references(scenario_path)
    used = reproduction["used"]
    predicted = [reproduction["predicted"][variant] for variant in used]
    panels = (
        ("gate_of_record", "Gate of record"),
        ("tavg_constant_dcp_half", "Dimer normalization only"),
        ("temp25_constant_dcp_dimer", "Temperature only (constant dCp)"),
        ("temp25_constant_dcp_half", "Temperature + normalization"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.2), sharex=True, sharey=True)
    for axis, (name, title) in zip(axes.flat, panels):
        observed = [references[variant][name] for variant in used]
        colors = ["#b44747" if variant.startswith("G93") else "#31688e" for variant in used]
        axis.scatter(observed, predicted, c=colors, s=48, edgecolor="white", linewidth=0.7)
        for variant, x_value, y_value in zip(used, observed, predicted):
            label_offset = {"G93S": (-30, 5), "G93A": (5, -13)}.get(variant, (4, 4))
            axis.annotate(variant, (x_value, y_value), xytext=label_offset,
                          textcoords="offset points", fontsize=8)
        axis.plot([-0.3, 7.3], [-0.3, 7.3], color="#555555", linestyle="--", linewidth=1)
        result = reproduction["metrics"][name]
        axis.set_title(f"{title}\nr={result['pearson']:.3f}; RMSE={result['rmse']:.3f}")
        axis.set_xlim(-0.3, 7.3)
        axis.set_ylim(-0.3, 7.3)
        axis.grid(alpha=0.2)
    fig.supxlabel("Experimental/reference delta-delta G (kcal/mol)")
    fig.supylabel("Frozen FEP prediction (kcal/mol)")
    fig.suptitle("Independent forensic reference-sensitivity reproduction")
    fig.tight_layout(rect=(0.02, 0.03, 1, 0.96))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", required=True, type=Path)
    parser.add_argument("--scenarios", type=Path,
                        default=Path("data/reference_sensitivity_scenarios.csv"))
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    result = reproduce(args.archive_root, args.scenarios)
    print(f"Used: {', '.join(result['used'])}; excluded: {', '.join(result['excluded'])}")
    print("| scenario | Pearson r | Spearman rho | RMSE | MUE |")
    print("|---|---:|---:|---:|---:|")
    for name, _label in SCENARIOS:
        row = result["metrics"][name]
        print(f"| {row['label']} | {row['pearson']:.3f} | {row['spearman']:.3f} | "
              f"{row['rmse']:.3f} | {row['mue']:.3f} |")
    if args.figure:
        plot(result, args.scenarios, args.figure)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
