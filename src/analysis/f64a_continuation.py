"""Fixed 54-fit CPU analysis of the complete, isolated 15 ns F64A continuation."""

from __future__ import annotations

import argparse
import json
from functools import partial
from pathlib import Path

import numpy as np
import yaml

from src.analysis import f64a_blocks as blocks
from src.analysis.f64a_sensitivity import fixed_trim_and_thin
from src.fep import f64a_continuation as pilot_module
from src.fep.analyze import decorrelate_window_with_diagnostics
from src.fep.f64a_extension import sha256_file


def load_design(path: Path, pilot: dict) -> dict:
    """Refuse unsafe output paths and incomplete or incompatible analysis plans."""
    design = yaml.safe_load(path.read_text())
    name = design["output_name"]
    if not isinstance(name, str) or Path(name).name != name or not name.endswith(".json"):
        raise ValueError("analysis output must be a JSON basename")
    tolerance = design["reproduction_atol"]
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("invalid reproduction tolerance")
    expected_policies = [
        {"label": "adaptive", "method": "adaptive_trim_and_thin"},
        {"label": "untrimmed", "method": "fixed_trim_and_thin", "trim_fraction": 0.},
        {"label": "fixed_25_percent", "method": "fixed_trim_and_thin", "trim_fraction": .25},
    ]
    if design["policies"] != expected_policies:
        raise ValueError("analysis policies differ from the registered three-policy plan")
    groups = [design["blocks"], design["cumulative"]]
    if [len(group) for group in groups] != [5, 3] or len(design["contrasts"]) != 3:
        raise ValueError("analysis must contain the registered 54-fit plan")
    for group in groups:
        if len({b["label"] for b in group}) != len(group):
            raise ValueError("duplicate analysis block")
        for block in group:
            start, stop = block["start"], block["stop"]
            if (type(start) is not int or type(stop) is not int
                    or not 0 <= start < stop <= pilot["target_samples_per_window"]):
                raise ValueError("invalid analysis slice")
    by_label = {b["label"]: b for b in design["blocks"]}
    if "block_6_9_ns" not in by_label or "cumulative_9_ns" not in {b["label"] for b in design["cumulative"]}:
        raise ValueError("missing v1 reproduction blocks")
    for contrast in design["contrasts"]:
        if (contrast["earlier"] not in by_label or contrast["later"] not in by_label
                or by_label[contrast["earlier"]]["stop"] > by_label[contrast["later"]]["start"]):
            raise ValueError("contrast must compare ordered, nonoverlapping blocks")
    return design


def check_reproduction(records: list[dict], reference: dict, design: dict) -> None:
    """Require adaptive 9 ns cumulative and 6–9 ns baseline metrics and counts."""
    for label, checkpoint in (("cumulative_9_ns", "9_ns"), ("block_6_9_ns", "final_3_ns")):
        for rep in range(3):
            record, = [r for r in records if r["policy"] == "adaptive" and r["block"] == label
                       and r["replicate"] == rep]
            expected, = [r for r in reference["records"] + reference["late_block_records"]
                         if r["checkpoint"] == checkpoint and r["replicate"] == rep]
            for field in blocks.REPRODUCTION_FIELDS:
                if not np.isclose(record[field], expected[field], rtol=0,
                                  atol=design["reproduction_atol"]):
                    raise ValueError(f"v1 reproduction failed: {label} r{rep} {field}")
            if record["n_raw"] != expected["n_raw"] or record["n_retained"] != expected["n_independent"]:
                raise ValueError("v1 reproduction failed: selected counts")


def analyze(config: str | Path = pilot_module.DEFAULT_CONFIG) -> Path:
    """Consume completed outputs only, preserving all policies and failed-fit visibility."""
    cfg, pilot = pilot_module.load_pilot(config)
    analysis_path = pilot_module.repo_path(pilot["analysis_config"])
    design = load_design(analysis_path, pilot)
    output = Path(pilot["_output_root"]) / "production/analysis" / design["output_name"]
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    fep, stage = pilot_module.validate_outputs(config)
    reference_path = pilot_module.repo_path(pilot["reference_primary"])
    reference = json.loads(reference_path.read_text())
    identity_design = {"_config_path": str(analysis_path), "_reference_path": str(reference_path),
                       "_extra_identity_paths": [str(Path(__file__).resolve()),
                        str(Path(pilot_module.__file__).resolve()),
                        str(pilot_module.repo_path("src/analysis/f64a_sensitivity.py")),
                        pilot["_config_path"], pilot["_base_config_path"]]}
    identity = blocks._analysis_identity(identity_design)
    paths = [fep / "F64A/folded" / f"w{w}_r{r}.npz" for w in range(20) for r in range(3)]
    paths += [fep.parent / "stage_manifest.json", Path(pilot["_output_root"]) / "admission.json", reference_path]
    paths += [fep.parent / "source/F64A/folded" / e["relative"] for e in stage["files"]]
    fingerprints = {str(p): sha256_file(p) for p in paths}
    records, summaries, contrasts, movements, local_changes = [], [], [], [], []
    for policy in design["policies"]:
        selector = (decorrelate_window_with_diagnostics if policy["method"] == "adaptive_trim_and_thin"
                    else partial(fixed_trim_and_thin, fraction=policy["trim_fraction"]))
        group = []
        for block in design["blocks"]:
            for rep in range(pilot["replicates"]):
                print(f"Solving {policy['label']} {block['label']} r{rep}", flush=True)
                result = blocks.solve_block(cfg, pilot, fep, block, rep, selector=selector)
                result.update(policy=policy["label"], kind="disjoint_block")
                group.append(result)
        a, b, c = blocks.summaries_and_contrasts(group, design)
        for source, target in ((a, summaries), (b, contrasts), (c, movements)):
            target.extend({"policy": policy["label"], **r} for r in source)
        for contrast in design["contrasts"]:
            for rep in range(3):
                earlier, = [r for r in group if r["block"] == contrast["earlier"] and r["replicate"] == rep]
                later, = [r for r in group if r["block"] == contrast["later"] and r["replicate"] == rep]
                for u, v in zip(earlier["adjacent_diagnostics"], later["adjacent_diagnostics"]):
                    local_changes.append({"policy": policy["label"], "contrast": contrast["label"],
                        "replicate": rep, "lower_window": u["lower_window"], "upper_window": u["upper_window"],
                        **{key + "_change": v[key] - u[key] for key in
                           ("forward_dg_kcal", "reverse_dg_kcal", "signed_discrepancy_kcal")}})
        records.extend(group)
    cumulative = []
    for block in design["cumulative"]:
        for rep in range(3):
            print(f"Solving adaptive {block['label']} r{rep}", flush=True)
            record = blocks.solve_block(cfg, pilot, fep, block, rep)
            record.update(policy="adaptive", kind="dependent_cumulative")
            cumulative.append(record)
    cumulative_summaries, _, _ = blocks.summaries_and_contrasts(
        cumulative, {"blocks": design["cumulative"], "contrasts": []})
    records.extend(cumulative)
    check_reproduction(records, reference, design)
    if any(sha256_file(p) != fingerprints[str(p)] for p in paths):
        raise ValueError("completed input changed during CPU analysis")
    if blocks._analysis_identity(identity_design) != identity:
        raise ValueError("analysis code or runtime changed during execution")
    pilot_module.write_json(output, {"schema_version": 1, "analysis_id": design["analysis_id"],
        "pilot_id": pilot["pilot_id"], "simulation_protocol": pilot_module.v1._extension_protocol(pilot, 6000.),
        "simulation_code_identity": stage["code_identity"], "source_code_identity": stage["source_code_identity"],
        "analysis_identity": identity, "input_sha256": fingerprints, "design": design,
        "scope": "exploratory F64A folded-only time dependence; no convergence pass or experimental gate",
        "v1_reproduction_passed": True, "records": records, "summaries": summaries,
        "cumulative_summaries": cumulative_summaries, "contrasts": contrasts,
        "contrast_summaries": movements, "local_changes": local_changes,
        "limitations": ["Policies and cumulative fits reuse data, not independent confirmations.",
            "Disjoint blocks can remain correlated through slow modes.",
            "Thinning counts/uncertainties are conditional and cannot cover unvisited conformations.",
            "Local diagnostics do not identify a structural mechanism or causal endpoint decomposition."]})
    print(f"Wrote {output}")
    print(json.dumps({"contrast_summaries": movements}, indent=2))
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=pilot_module.DEFAULT_CONFIG)
    analyze(parser.parse_args().config)


if __name__ == "__main__":
    main()
