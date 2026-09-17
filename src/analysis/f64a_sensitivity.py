"""CPU-only selection-policy sensitivity for the completed F64A folded pilot."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from functools import partial
from pathlib import Path

import numpy as np

from src.analysis import f64a_blocks as blocks
from src.fep.analyze import decorrelate_window_with_diagnostics
from src.fep.f64a_extension import ROOT, sha256_file, validate_production_outputs

DEFAULT_CONFIG = ROOT / "pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml"


def fixed_trim_and_thin(u: np.ndarray, window: int, fraction: float) -> tuple[np.ndarray, dict]:
    """Apply a fixed cutoff; estimate both correlation times on the remaining suffix.

    Never invoke adaptive equilibration detection. Correlation thinning does not
    establish stationarity or independent conformational coverage.
    """
    if (u.ndim != 2 or u.shape[0] < 2 or not 0 <= window < u.shape[0]
            or u.shape[1] == 0 or not np.isfinite(u).all()):
        raise ValueError("invalid or non-finite sensitivity window")
    if not np.isfinite(fraction) or not 0 <= fraction < 1:
        raise ValueError("invalid fixed trim fraction")
    from pymbar import timeseries

    cutoff = int(np.floor(fraction * u.shape[1]))
    neighbor = window + 1 if window + 1 < u.shape[0] else window - 1
    own = u[window, cutoff:]
    difference = u[neighbor, cutoff:] - own
    if len(own) < 10:
        raise ValueError("fixed suffix too short for correlation estimation")
    g_self = float(timeseries.statistical_inefficiency(own))
    g_neighbor = float(timeseries.statistical_inefficiency(difference))
    g = max(g_self, g_neighbor)
    if not np.isfinite(g) or g < 1:
        raise ValueError("invalid correlation estimate")
    keep = timeseries.subsample_correlated_data(own, g=g)
    if not len(keep):
        raise ValueError("fixed trimming/thinning retained no samples")
    return u[:, cutoff:][:, keep], {
        "method": "fixed_trim_and_thin", "trim_fraction": fraction,
        "neighbor_state": neighbor, "n_raw": int(u.shape[1]),
        "t0_self": None, "t0_neighbor": None, "trim_start_column": cutoff,
        "g_self": g_self, "g_neighbor": g_neighbor, "subsample_g": g,
        "n_after_trim": int(len(own)), "n_retained": int(len(keep)),
        "first_selected_column": int(cutoff + keep[0]),
        "last_selected_column": int(cutoff + keep[-1]),
    }


def load_design(config: str | Path = DEFAULT_CONFIG) -> tuple[dict, dict, dict]:
    """Load the fixed three-policy comparison and register its frozen reference."""
    cfg, pilot, design = blocks.load_design(config)
    expected = [
        {"label": "adaptive", "method": "adaptive_trim_and_thin"},
        {"label": "untrimmed", "method": "fixed_trim_and_thin", "trim_fraction": 0.0},
        {"label": "fixed_25_percent", "method": "fixed_trim_and_thin", "trim_fraction": 0.25},
    ]
    if design["policies"] != expected:
        raise ValueError("sensitivity requires the fixed adaptive/0%/25% policies")
    if design["output_name"] != "f64a_selection_sensitivity_v1.json":
        raise ValueError("sensitivity requires its isolated output name")
    if [(b["label"], b["start"], b["stop"]) for b in design["blocks"]] != [
        ("block_6_7p5_ns", 6001, 7501), ("block_7p5_9_ns", 7501, 9001),
    ]:
        raise ValueError("sensitivity requires the two registered final halves")
    if design["contrasts"] != [{"label": "final_half_to_half",
                                "earlier": "block_6_7p5_ns", "later": "block_7p5_9_ns"}]:
        raise ValueError("sensitivity requires the registered paired contrast")
    design["_block_reference_path"] = str(blocks._repo_path(design["block_reference"]))
    design["_extra_identity_paths"] = [str(Path(__file__).resolve()),
                                       design["_block_reference_path"]]
    return cfg, pilot, design


def check_adaptive_reproduction(records: list[dict], reference: dict, design: dict) -> None:
    """Reproduce all six adaptive fits, including exact per-window sample selection."""
    for record in records:
        if record["policy"] != "adaptive":
            continue
        matches = [r for r in reference["records"] if r["block"] == record["block"]
                   and r["replicate"] == record["replicate"]]
        if len(matches) != 1:
            raise ValueError("missing or ambiguous adaptive reference")
        expected = matches[0]
        for field in blocks.REPRODUCTION_FIELDS:
            if not np.isclose(record[field], expected[field], rtol=0,
                              atol=design["reproduction_atol"]):
                raise ValueError(f"adaptive reproduction failed: {field}")
        for field in ("n_raw", "n_retained", "window_diagnostics"):
            if record[field] != expected[field]:
                raise ValueError(f"adaptive reproduction failed: {field}")


def analyze_sensitivity(config: str | Path = DEFAULT_CONFIG) -> Path:
    """Validate immutable completed inputs and write an isolated 18-fit CPU report."""
    cfg, pilot, design = load_design(config)
    output = Path(pilot["_output_root"]) / "production/analysis" / design["output_name"]
    if output.exists():
        raise FileExistsError(f"refusing to overwrite sensitivity analysis {output}")
    fep, stage = validate_production_outputs(design["_pilot_config_path"])
    primary = fep.parent / "analysis/f64a_extension.json"
    block_report = fep.parent / "analysis/f64a_blocks_v1.json"
    frozen_primary = Path(design["_reference_path"])
    frozen_blocks = Path(design["_block_reference_path"])
    reference = json.loads(frozen_primary.read_text())
    block_reference = json.loads(frozen_blocks.read_text())
    if json.loads(primary.read_text()) != reference:
        raise ValueError("cluster primary differs from the preserved report")
    if json.loads(block_report.read_text()) != block_reference:
        raise ValueError("cluster blocks differ from the preserved report")
    if (stage["code_identity"] != reference["code_identity"]
            or block_reference["simulation_code_identity"] != stage["code_identity"]
            or block_reference["simulation_protocol"] != reference["protocol"]
            or block_reference["reproduction_passed"] is not True):
        raise ValueError("incompatible simulation/block reference identities")
    identity = blocks._analysis_identity(design)
    paths = [fep / "F64A/folded" / f"w{w}_r{r}.npz"
             for w in range(pilot["lambda_windows"]) for r in range(pilot["replicates"])]
    paths += [fep.parent / "stage_manifest.json", primary, block_report,
              frozen_primary, frozen_blocks]
    fingerprints = {str(p): sha256_file(p) for p in paths}
    previous_npz = {Path(p).name: h for p, h in block_reference["input_sha256"].items()
                    if p.endswith(".npz")}
    for path in paths:
        if path.suffix == ".npz" and previous_npz.get(path.name) != fingerprints[str(path)]:
            raise ValueError(f"production NPZ differs from the frozen block input: {path.name}")
    records, summaries, contrasts, contrast_summaries, local_changes = [], [], [], [], []
    for policy in design["policies"]:
        selector = (decorrelate_window_with_diagnostics if policy["label"] == "adaptive"
                    else partial(fixed_trim_and_thin, fraction=policy["trim_fraction"]))
        group = []
        for block in design["blocks"]:
            for rep in range(pilot["replicates"]):
                print(f"Solving {policy['label']} {block['label']} r{rep}", flush=True)
                record = blocks.solve_block(cfg, pilot, fep, block, rep, selector=selector)
                record["policy"] = policy["label"]
                group.append(record)
        a, b, c = blocks.summaries_and_contrasts(group, design)
        for source, target in ((a, summaries), (b, contrasts), (c, contrast_summaries)):
            target.extend({"policy": policy["label"], **r} for r in source)
        for rep in range(pilot["replicates"]):
            earlier, later = [r for r in group if r["replicate"] == rep]
            for u, v in zip(earlier["adjacent_diagnostics"], later["adjacent_diagnostics"]):
                local_changes.append({
                    "policy": policy["label"], "replicate": rep,
                    "lower_window": u["lower_window"], "upper_window": u["upper_window"],
                    "earlier_signed_discrepancy_kcal": u["signed_discrepancy_kcal"],
                    "later_signed_discrepancy_kcal": v["signed_discrepancy_kcal"],
                    "signed_discrepancy_change_kcal":
                        v["signed_discrepancy_kcal"] - u["signed_discrepancy_kcal"],
                    "forward_dg_change_kcal": v["forward_dg_kcal"] - u["forward_dg_kcal"],
                    "reverse_dg_change_kcal": v["reverse_dg_kcal"] - u["reverse_dg_kcal"],
                })
        records.extend(group)
    check_adaptive_reproduction(records, block_reference, design)
    if any(sha256_file(p) != fingerprints[str(p)] for p in paths):
        raise ValueError("completed input changed during sensitivity analysis")
    if blocks._analysis_identity(design) != identity:
        raise ValueError("analysis code or runtime changed during execution")
    report = {
        "schema_version": 1, "analysis_id": design["analysis_id"],
        "pilot_id": reference["pilot_id"], "simulation_protocol": reference["protocol"],
        "simulation_code_identity": stage["code_identity"], "analysis_identity": identity,
        "input_sha256": fingerprints,
        "design": {k: v for k, v in design.items() if not k.startswith("_")},
        "scope": "exploratory folded-only CPU sensitivity; no new sampling or gate",
        "adaptive_reproduction_passed": True, "records": records, "summaries": summaries,
        "contrasts": contrasts, "contrast_summaries": contrast_summaries,
        "local_hysteresis_changes": local_changes,
        "limitations": [
            "Policies share raw data; differences are not independent replications.",
            "Fixed 25% trimming is an exploratory sensitivity, not an equilibration criterion.",
            "Correlation estimates from nonstationary intervals can be unreliable.",
            "Conditional uncertainties cannot cover unsampled conformations or model/reference bias.",
            "Local hysteresis changes do not decompose endpoint drift or identify a structural mode.",
        ],
    }
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=output.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(serialized)
    try:
        os.link(temporary, output)  # Atomic publication, refusing even a concurrent overwrite.
    finally:
        temporary.unlink()
    print(f"Wrote {output}")
    print(json.dumps({"contrast_summaries": contrast_summaries}, indent=2))
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    analyze_sensitivity(parser.parse_args().config)


if __name__ == "__main__":
    main()
