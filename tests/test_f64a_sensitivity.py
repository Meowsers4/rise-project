"""Synthetic contracts, not new scientific F64A results."""

import copy
import json
import shutil
import subprocess

import numpy as np
import pytest
import yaml

from src.analysis import f64a_sensitivity as sensitivity
from src.fep.f64a_extension import ROOT, sha256_file


@pytest.mark.parametrize("fraction,cutoff", [(0., 0), (.25, 10)])
def test_fixed_policy_never_detects_equilibration_and_reestimates_suffix_g(
        monkeypatch, fraction, cutoff):
    from pymbar import timeseries

    u = np.vstack([np.arange(40.), np.arange(40.) * 3])
    seen = []

    def inefficiency(series):
        seen.append(series.copy())
        return 2. if len(seen) == 1 else 4.

    monkeypatch.setattr(timeseries, "detect_equilibration",
                        lambda *a: pytest.fail("fixed policy cannot auto-trim"))
    monkeypatch.setattr(timeseries, "statistical_inefficiency", inefficiency)
    monkeypatch.setattr(timeseries, "subsample_correlated_data",
                        lambda series, g: [0, 4, 8] if g == 4. else [])
    selected, diagnostic = sensitivity.fixed_trim_and_thin(u, 0, fraction)
    assert np.array_equal(seen[0], u[0, cutoff:])
    assert np.array_equal(seen[1], (u[1] - u[0])[cutoff:])
    assert np.array_equal(selected, u[:, [cutoff, cutoff + 4, cutoff + 8]])
    assert diagnostic["trim_start_column"] == cutoff
    assert diagnostic["t0_self"] is None
    assert diagnostic["last_selected_column"] == cutoff + 8


@pytest.mark.parametrize("fraction", [-.1, 1., float("nan")])
def test_invalid_fixed_trim_refused(fraction):
    with pytest.raises(ValueError, match="fraction"):
        sensitivity.fixed_trim_and_thin(np.ones((2, 40)), 0, fraction)


def test_failed_correlation_is_not_silently_replaced(monkeypatch):
    from pymbar import timeseries

    monkeypatch.setattr(timeseries, "statistical_inefficiency", lambda *a: float("nan"))
    with pytest.raises(ValueError, match="correlation"):
        sensitivity.fixed_trim_and_thin(np.ones((2, 40)), 0, .25)


def test_registered_design_and_frozen_block_reference():
    _, _, design = sensitivity.load_design()
    assert len(design["policies"]) == 3
    reference = json.loads(open(design["_block_reference_path"]).read())
    assert reference["analysis_identity"]["git_commit"].startswith("0ba895ab")
    assert sum(len(r["window_diagnostics"]) for r in reference["records"]) == 300
    records = [{**r, "policy": "adaptive"} for r in reference["records"]
               if r["block"] in {b["label"] for b in design["blocks"]}]
    assert len(records) == 6
    sensitivity.check_adaptive_reproduction(records, reference, design)
    records[0] = copy.deepcopy(records[0])
    records[0]["window_diagnostics"][0]["n_retained"] -= 1
    with pytest.raises(ValueError, match="window_diagnostics"):
        sensitivity.check_adaptive_reproduction(records, reference, design)


def _mock_pipeline(tmp_path, monkeypatch):
    cfg, pilot, design = sensitivity.load_design()
    pilot = {**pilot, "_output_root": str(tmp_path)}
    primary = json.loads(open(design["_reference_path"]).read())
    reference = json.loads(open(design["_block_reference_path"]).read())
    fep = tmp_path / "production/fep"
    folder = fep / "F64A/folded"
    folder.mkdir(parents=True)
    for w in range(20):
        for r in range(3):
            (folder / f"w{w}_r{r}.npz").write_bytes(b"mock validated input")
    reference["input_sha256"] = {str(p): sha256_file(p) for p in folder.iterdir()}
    frozen_blocks = tmp_path / "frozen_blocks.json"
    frozen_blocks.write_text(json.dumps(reference))
    design["_block_reference_path"] = str(frozen_blocks)
    analysis = tmp_path / "production/analysis"
    analysis.mkdir()
    (analysis / "f64a_extension.json").write_text(json.dumps(primary))
    (analysis / "f64a_blocks_v1.json").write_text(json.dumps(reference))
    (tmp_path / "production/stage_manifest.json").write_text("{}")
    monkeypatch.setattr(sensitivity, "load_design", lambda config: (cfg, pilot, design))
    monkeypatch.setattr(sensitivity, "validate_production_outputs",
                        lambda config: (fep, {"code_identity": primary["code_identity"]}))
    monkeypatch.setattr(sensitivity.blocks, "_analysis_identity", lambda design: {"test": True})

    def solve(cfg, pilot, fep, block, rep, selector):
        record = copy.deepcopy(next(r for r in reference["records"]
                                    if r["block"] == block["label"] and r["replicate"] == rep))
        for pair in record["adjacent_diagnostics"]:
            pair["forward_dg_kcal"] = 0.
            pair["reverse_dg_kcal"] = -pair["signed_discrepancy_kcal"]
        return record

    monkeypatch.setattr(sensitivity.blocks, "solve_block", solve)
    return fep, analysis


def test_pipeline_18_fits_preserves_previous_reports_and_refuses_overwrite(tmp_path, monkeypatch):
    _, analysis = _mock_pipeline(tmp_path, monkeypatch)
    prior = {p: p.read_bytes() for p in analysis.iterdir()}
    result = json.loads(sensitivity.analyze_sensitivity().read_text())
    assert len(result["records"]) == 18
    assert len(result["contrast_summaries"]) == 3
    assert len(result["local_hysteresis_changes"]) == 171
    assert result["adaptive_reproduction_passed"]
    assert all(p.read_bytes() == content for p, content in prior.items())
    with pytest.raises(FileExistsError):
        sensitivity.analyze_sensitivity()


@pytest.mark.parametrize("failure", ["reference", "fit", "mutation", "reproduction", "old_hash"])
def test_failure_writes_no_partial_report(tmp_path, monkeypatch, failure):
    fep, analysis = _mock_pipeline(tmp_path, monkeypatch)
    original = sensitivity.blocks.solve_block
    if failure == "reference":
        (analysis / "f64a_blocks_v1.json").write_text("{}")
    elif failure == "old_hash":
        (fep / "F64A/folded/w0_r0.npz").write_bytes(b"changed before run")
    else:
        def solve(*args, **kwargs):
            if failure == "fit":
                raise ValueError("failed fit")
            record = original(*args, **kwargs)
            if failure == "mutation":
                (fep / "F64A/folded/w0_r0.npz").write_bytes(b"changed input")
            if failure == "reproduction":
                record["dg_kcal"] += .1
            return record
        monkeypatch.setattr(sensitivity.blocks, "solve_block", solve)
    with pytest.raises(ValueError):
        sensitivity.analyze_sensitivity()
    assert not (analysis / "f64a_selection_sensitivity_v1.json").exists()


def test_snakemake_sensitivity_cpu_only_and_missing_report_fatal(tmp_path):
    executable = shutil.which("snakemake")
    if not executable:
        pytest.skip("Snakemake unavailable")
    for directory in ("config", "workflow", "pilots/f64a_sampling_extension_v1"):
        (tmp_path / directory).mkdir(parents=True)
    for path in ("config/pipeline.yaml", "workflow/Snakefile",
                 "pilots/f64a_sampling_extension_v1/cpu_blocks.yaml",
                 "pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml"):
        (tmp_path / path).write_bytes((ROOT / path).read_bytes())
    pilot = yaml.safe_load((ROOT / "pilots/f64a_sampling_extension_v1/config.yaml").read_text())
    pilot["output_root"] = "mock_completed"
    (tmp_path / "pilots/f64a_sampling_extension_v1/config.yaml").write_text(yaml.safe_dump(pilot))
    for directory in ("src", "docs", "data"):
        (tmp_path / directory).symlink_to(ROOT / directory, target_is_directory=True)
    completed = tmp_path / "mock_completed/production"
    (completed / "analysis").mkdir(parents=True)
    (completed / "stage_manifest.json").write_text("{}")
    for report in ("f64a_extension.json", "f64a_blocks_v1.json"):
        (completed / "analysis" / report).write_text("{}")
    folder = completed / "fep/F64A/folded"
    folder.mkdir(parents=True)
    for w in range(20):
        for r in range(3):
            (folder / f"w{w}_r{r}.npz").write_bytes(b"mock completed input")
    command = [executable, "f64a_extension_sensitivity", "-s", "workflow/Snakefile", "-n",
               "--cores", "1", "--runtime-source-cache-path", str(tmp_path / "cache"),
               "--allowed-rules", "f64a_extension_sensitivity"]
    run = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    output = run.stdout + run.stderr
    assert run.returncode == 0, output
    assert "rule f64a_extension_sensitivity:" in output
    assert "rule f64a_extension_window:" not in output
    assert "rule f64a_extension_blocks:" not in output
    (completed / "analysis/f64a_blocks_v1.json").unlink()
    run = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    assert run.returncode != 0
    assert "will not schedule GPU work" in run.stdout + run.stderr
