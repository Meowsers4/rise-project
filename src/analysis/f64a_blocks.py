"""CPU-only exploratory disjoint-block diagnostic for the completed F64A pilot."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path

import numpy as np
import yaml

from src.fep.analyze import (
    _KB_KCAL,
    decorrelate_window_with_diagnostics,
    solve_leg_mbar,
)
from src.fep.f64a_extension import (
    ROOT,
    adjacent_discrepancies,
    load_pilot,
    sha256_file,
    validate_production_outputs,
)

DEFAULT_CONFIG = ROOT / "pilots/f64a_sampling_extension_v1/cpu_blocks.yaml"
REPRODUCTION_FIELDS = (
    "dg_kcal", "mbar_error_kcal", "net_hysteresis_kcal",
    "sum_absolute_hysteresis_kcal", "max_absolute_hysteresis_kcal",
    "min_adjacent_overlap",
)


def _repo_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def load_design(config: str | Path = DEFAULT_CONFIG) -> tuple[dict, dict, dict]:
    """Validate half-open intervals and reject overlapping paired contrasts."""
    config = Path(config).resolve()
    design = yaml.safe_load(config.read_text())
    cfg, pilot = load_pilot(_repo_path(design["pilot_config"]))
    if design["method"] != "adaptive_trim_and_thin":
        raise ValueError("unsupported block-selection method")
    tolerance = float(design["reproduction_atol"])
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("invalid numerical reproduction tolerance")
    output_name = design["output_name"]
    if Path(output_name).name != output_name or not output_name.endswith(".json"):
        raise ValueError("output must be a JSON basename inside the pilot analysis directory")
    if output_name == "f64a_extension.json":
        raise ValueError("refusing to replace the primary pilot analysis")
    blocks = {}
    for block in design["blocks"]:
        start, stop = block["start"], block["stop"]
        if (type(start) is not int or type(stop) is not int
                or not 0 <= start < stop <= pilot["target_samples_per_window"]):
            raise ValueError(f"invalid block slice: {block}")
        if block["label"] in blocks:
            raise ValueError("duplicate block label")
        blocks[block["label"]] = block
    if not blocks:
        raise ValueError("empty block design")
    if len({c["label"] for c in design["contrasts"]}) != len(design["contrasts"]):
        raise ValueError("duplicate contrast label")
    for contrast in design["contrasts"]:
        earlier, later = blocks[contrast["earlier"]], blocks[contrast["later"]]
        if earlier["stop"] > later["start"]:
            raise ValueError(f"contrast is overlapping or reversed: {contrast}")
    design["_config_path"] = str(config)
    design["_pilot_config_path"] = str(_repo_path(design["pilot_config"]))
    design["_reference_path"] = str(_repo_path(design["reference_report"]))
    return cfg, pilot, design


def _analysis_identity(design: dict) -> dict:
    paths = [Path(__file__).resolve(), ROOT / "src/fep/analyze.py",
             ROOT / "src/fep/f64a_extension.py", Path(design["_config_path"]),
             Path(design["_reference_path"])]
    relative_paths = [str(path.relative_to(ROOT)) for path in paths]
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                            capture_output=True, text=True).stdout.strip()
    status = subprocess.run(["git", "status", "--porcelain", "--", *relative_paths],
                            cwd=ROOT, check=True, capture_output=True,
                            text=True).stdout.splitlines()
    if status:
        raise RuntimeError(f"commit the CPU analyzer/config/reference before execution: {status}")
    versions = {}
    for package in ("numpy", "pymbar", "jax", "jaxlib", "PyYAML"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return {
        "git_commit": commit, "analysis_paths_clean": True,
        "source_sha256": {name: sha256_file(path)
                          for name, path in zip(relative_paths, paths)},
        "python": platform.python_version(), "packages": versions,
    }


def solve_block(cfg: dict, pilot: dict, fep_dir: Path, block: dict, rep: int) -> dict:
    """Solve one block/repeat and retain all per-window and adjacent diagnostics."""
    start, stop = block["start"], block["stop"]
    windows, diagnostics = [], []
    stride = float(pilot["energy_stride_ps"])
    baseline_start = float(pilot["baseline_total_time_ps"]) - pilot["baseline_ns"] * 1000
    for window in range(pilot["lambda_windows"]):
        with np.load(fep_dir / "F64A/folded" / f"w{window}_r{rep}.npz") as data:
            u = np.asarray(data["u_kn_window"])
        if u.shape != (pilot["lambda_windows"], pilot["target_samples_per_window"]):
            raise ValueError("short or incompatible block input; refusing silent slicing")
        selected, choices = decorrelate_window_with_diagnostics(u[:, start:stop], window)
        windows.append(selected)
        record = {"window": window, **choices}
        for name in ("trim_start_column", "first_selected_column", "last_selected_column"):
            raw_column = start + choices[name]
            record[f"{name}_global"] = raw_column
            record[f"{name}_retained_time_ps"] = raw_column * stride
            record[f"{name}_simulation_time_ps"] = baseline_start + raw_column * stride
        diagnostics.append(record)
    n_k = np.asarray([u.shape[1] for u in windows])
    solved = solve_leg_mbar(np.concatenate(windows, axis=1), n_k)
    kt = _KB_KCAL * cfg["fep"]["temperature_K"]
    local = adjacent_discrepancies(windows, kt)
    overlap = np.asarray(solved["overlap"])
    return {
        "replicate": rep, "block": block["label"],
        "raw_column_start": start, "raw_column_stop": stop,
        "n_raw": int((stop - start) * pilot["lambda_windows"]),
        "n_retained": int(n_k.sum()),
        "dg_kcal": solved["dg_kT"] * kt,
        "mbar_error_kcal": solved["ddg_kT"] * kt,
        "net_hysteresis_kcal": float(abs(local.sum())),
        "sum_absolute_hysteresis_kcal": float(abs(local).sum()),
        "max_absolute_hysteresis_kcal": float(abs(local).max()),
        "min_adjacent_overlap": solved["min_adjacent"],
        "overlap_matrix": solved["overlap"], "solver_notes": solved["solver_notes"],
        "window_diagnostics": diagnostics,
        "adjacent_diagnostics": [
            {"lower_window": k, "upper_window": k + 1,
             "signed_discrepancy_kcal": float(value),
             "overlap_lower_to_upper": float(overlap[k, k + 1]),
             "overlap_upper_to_lower": float(overlap[k + 1, k])}
            for k, value in enumerate(local)
        ],
    }


def check_reproduction(records: list[dict], reference: dict, design: dict) -> None:
    """Require the two already-reported blocks to reproduce, without re-gating."""
    original = reference["records"] + reference["late_block_records"]
    for block in design["blocks"]:
        if "reference" not in block:
            continue
        for record in [r for r in records if r["block"] == block["label"]]:
            match = [r for r in original if r["checkpoint"] == block["reference"]
                     and r["replicate"] == record["replicate"]]
            if len(match) != 1:
                raise ValueError("missing or ambiguous reference block")
            expected = match[0]
            for field in REPRODUCTION_FIELDS:
                if not np.isclose(record[field], expected[field], rtol=0,
                                  atol=design["reproduction_atol"]):
                    raise ValueError(f"reference reproduction failed: {block['label']} "
                                     f"r{record['replicate']} {field}")
            if (record["n_raw"] != expected["n_raw"]
                    or record["n_retained"] != expected["n_independent"]):
                raise ValueError("reference reproduction failed: sample counts")


def summaries_and_contrasts(records: list[dict], design: dict) -> tuple[list, list, list]:
    """Summarize three repeat estimates and paired movements descriptively."""
    summaries, contrasts, contrast_summaries = [], [], []
    for block in design["blocks"]:
        group = [r for r in records if r["block"] == block["label"]]
        values = np.asarray([r["dg_kcal"] for r in group])
        sem = float(values.std(ddof=1) / np.sqrt(len(values)))
        propagated = float(np.sqrt(sum(r["mbar_error_kcal"] ** 2 for r in group)) / len(group))
        summaries.append({
            "block": block["label"], "mean_dg_kcal": float(values.mean()),
            "replicate_range_kcal": float(np.ptp(values)), "replicate_sem_kcal": sem,
            "propagated_mbar_error_kcal": propagated, "uncertainty_kcal": max(sem, propagated),
        })
    for contrast in design["contrasts"]:
        earlier = {r["replicate"]: r for r in records if r["block"] == contrast["earlier"]}
        later = {r["replicate"]: r for r in records if r["block"] == contrast["later"]}
        changes = []
        for rep in sorted(earlier):
            delta = later[rep]["dg_kcal"] - earlier[rep]["dg_kcal"]
            changes.append(delta)
            contrasts.append({
                **contrast, "replicate": rep, "signed_movement_kcal": delta,
                "absolute_movement_kcal": abs(delta), "raw_columns_disjoint": True,
                "earlier_mbar_error_kcal": earlier[rep]["mbar_error_kcal"],
                "later_mbar_error_kcal": later[rep]["mbar_error_kcal"],
            })
        contrast_summaries.append({
            **contrast, "mean_signed_movement_kcal": float(np.mean(changes)),
            "paired_movement_sem_kcal": float(np.std(changes, ddof=1) / np.sqrt(len(changes))),
            "max_absolute_movement_kcal": float(np.max(np.abs(changes))),
        })
    return summaries, contrasts, contrast_summaries


def analyze_blocks(config: str | Path = DEFAULT_CONFIG) -> Path:
    """Validate immutable completed inputs, perform CPU fits, and write a new report."""
    cfg, pilot, design = load_design(config)
    output = Path(pilot["_output_root"]) / "production/analysis" / design["output_name"]
    if output.exists():
        raise FileExistsError(f"refusing to overwrite block analysis {output}")
    fep_dir, stage = validate_production_outputs(design["_pilot_config_path"])
    primary_path = fep_dir.parent / "analysis/f64a_extension.json"
    reference_path = Path(design["_reference_path"])
    reference = json.loads(reference_path.read_text())
    if json.loads(primary_path.read_text()) != reference:
        raise ValueError("cluster primary report differs from the preserved supplied report")
    if stage["code_identity"] != reference["code_identity"]:
        raise ValueError("primary report and staged simulation identities differ")
    identity = _analysis_identity(design)
    paths = [fep_dir / "F64A/folded" / f"w{w}_r{r}.npz"
             for w in range(pilot["lambda_windows"]) for r in range(pilot["replicates"])]
    paths += [fep_dir.parent / "stage_manifest.json", primary_path, reference_path]
    fingerprints = {str(path): sha256_file(path) for path in paths}
    records = []
    for block in design["blocks"]:
        for rep in range(pilot["replicates"]):
            print(f"Solving {block['label']} r{rep}", flush=True)
            records.append(solve_block(cfg, pilot, fep_dir, block, rep))
    check_reproduction(records, reference, design)
    summaries, contrasts, contrast_summaries = summaries_and_contrasts(records, design)
    if any(sha256_file(path) != fingerprints[str(path)] for path in paths):
        raise ValueError("completed input changed during CPU analysis")
    if _analysis_identity(design) != identity:
        raise ValueError("analysis code or runtime identity changed during execution")
    report = {
        "schema_version": 1, "analysis_id": design["analysis_id"],
        "pilot_id": reference["pilot_id"], "simulation_protocol": reference["protocol"],
        "pilot_config_sha256": reference["pilot_config_sha256"],
        "base_config_sha256": reference["base_config_sha256"],
        "scope": "exploratory F64A folded-only; no convergence gate or experimental comparison",
        "simulation_code_identity": stage["code_identity"],
        "analysis_identity": identity, "input_sha256": fingerprints,
        "design": {k: v for k, v in design.items() if not k.startswith("_")},
        "reproduction_passed": True, "records": records, "summaries": summaries,
        "contrasts": contrasts, "contrast_summaries": contrast_summaries,
        "limitations": [
            "Disjoint time blocks can remain correlated through slow unsampled modes.",
            "Adaptive trimming changes between blocks; retained counts are heuristic.",
            "The final block and its two halves overlap; do not treat all five blocks as independent.",
            "Local hysteresis is not a decomposition of endpoint drift or a structural mechanism.",
        ],
    }
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(serialized)
    temporary.replace(output)
    print(f"Wrote {output}")
    print(json.dumps({"summaries": summaries, "contrast_summaries": contrast_summaries}, indent=2))
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    analyze_blocks(args.config)


if __name__ == "__main__":
    main()
