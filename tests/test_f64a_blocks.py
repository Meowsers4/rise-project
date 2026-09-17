"""Contracts for the exploratory CPU follow-up; synthetic fixtures are not results."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess

import numpy as np
import pytest
import yaml

import src.analysis.f64a_blocks as blocks
from src.fep.analyze import decorrelate_window, decorrelate_window_with_diagnostics
from src.fep.f64a_extension import ROOT, sha256_file


def _reference():
    return json.loads((ROOT / "docs/f64a_extension_review/f64a_extension.json").read_text())


def _design_file(tmp_path, change):
    design = yaml.safe_load(blocks.DEFAULT_CONFIG.read_text())
    change(design)
    path = tmp_path / "blocks.yaml"
    path.write_text(yaml.safe_dump(design))
    return path


def test_disjoint_blocks_partition_without_duplicating_boundary_columns():
    _, pilot, design = blocks.load_design()
    first, middle, final, half_a, half_b = design["blocks"]
    assert [b["stop"] - b["start"] for b in design["blocks"]] == [3001, 3000, 3000, 1500, 1500]
    assert first["start"] == 0
    assert first["stop"] == middle["start"]
    assert middle["stop"] == final["start"]
    assert final["stop"] == pilot["target_samples_per_window"]
    assert half_a["start"] == final["start"]
    assert half_a["stop"] == half_b["start"]
    assert half_b["stop"] == final["stop"]


@pytest.mark.parametrize("change, error", [
    (lambda d: d["blocks"][1].update(start=3000), "overlapping"),
    (lambda d: d["blocks"][0].update(stop=9002), "invalid block"),
    (lambda d: d["blocks"][0].update(start=0.0), "invalid block"),
    (lambda d: d.update(output_name="../primary.json"), "basename"),
    (lambda d: d.update(output_name="f64a_extension.json"), "primary"),
])
def test_design_rejects_unsafe_or_ambiguous_slices(tmp_path, change, error):
    with pytest.raises(ValueError, match=error):
        blocks.load_design(_design_file(tmp_path, change))


def test_diagnostic_selection_preserves_historical_algorithm(monkeypatch):
    from pymbar import timeseries

    u = np.vstack([np.arange(30.), np.arange(30.) * 3])

    def estimate(series):
        return (2, 3.0, 9.0) if np.array_equal(series, u[0]) else (5, 7.0, 4.0)

    def subsample(series, g):
        assert g == 7.0
        assert np.array_equal(series, u[0, 5:])
        return [0, 7, 14]

    monkeypatch.setattr(timeseries, "detect_equilibration", estimate)
    monkeypatch.setattr(timeseries, "subsample_correlated_data", subsample)
    selected, diagnostics = decorrelate_window_with_diagnostics(u, 0)
    expected = u[:, 5:][:, [0, 7, 14]]
    assert np.array_equal(selected, expected)
    assert np.array_equal(decorrelate_window(u, 0), expected)
    assert diagnostics["t0_self"] == 2
    assert diagnostics["t0_neighbor"] == 5
    assert diagnostics["trim_start_column"] == 5
    assert diagnostics["n_retained"] == 3
    assert diagnostics["last_selected_column"] == 19


def test_real_timeseries_selection_unchanged():
    from pymbar import timeseries

    rng = np.random.default_rng(20260916)
    u = rng.normal(size=(3, 120))
    u_self, du = u[1], u[2] - u[1]
    t0, g_self, _ = timeseries.detect_equilibration(u_self)
    t0_du, g_du, _ = timeseries.detect_equilibration(du)
    t0, g = max(t0, t0_du), max(g_self, g_du)
    keep = timeseries.subsample_correlated_data(u_self[t0:], g=g)
    selected, diagnostics = decorrelate_window_with_diagnostics(u, 1)
    assert np.array_equal(selected, u[:, t0:][:, keep])
    assert diagnostics["n_retained"] == len(keep)


def test_short_block_and_nonfinite_guard():
    u = np.zeros((2, 9))
    selected, diagnostic = decorrelate_window_with_diagnostics(u, 0)
    assert np.array_equal(selected, u)
    assert diagnostic["method"] == "short_block_no_thinning"
    u[1, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        decorrelate_window_with_diagnostics(u, 0)


def test_solver_records_global_times_and_all_adjacent_pairs(tmp_path, monkeypatch):
    cfg, pilot, _ = blocks.load_design()
    pilot = {**pilot, "lambda_windows": 2, "target_samples_per_window": 12}
    folder = tmp_path / "F64A/folded"
    folder.mkdir(parents=True)
    for window in range(2):
        np.savez(folder / f"w{window}_r0.npz", u_kn_window=np.zeros((2, 12)))
    monkeypatch.setattr(blocks, "solve_leg_mbar", lambda *args: {
        "dg_kT": 1., "ddg_kT": .1, "overlap": [[.9, .1], [.2, .8]],
        "min_adjacent": .1, "solver_notes": [],
    })
    result = blocks.solve_block(cfg, pilot, tmp_path, {"label": "tiny", "start": 4, "stop": 9}, 0)
    assert len(result["window_diagnostics"]) == 2
    assert len(result["adjacent_diagnostics"]) == 1
    diagnostic = result["window_diagnostics"][0]
    assert diagnostic["first_selected_column_global"] == 4
    assert diagnostic["first_selected_column_retained_time_ps"] == 4
    assert diagnostic["first_selected_column_simulation_time_ps"] == 504
    assert result["n_retained"] == 10
    np.savez(folder / "w0_r0.npz", u_kn_window=np.zeros((2, 8)))
    with pytest.raises(ValueError, match="silent slicing"):
        blocks.solve_block(cfg, pilot, tmp_path, {"label": "tiny", "start": 4, "stop": 9}, 0)


def _records_for_design(design, reference):
    records = []
    for block in design["blocks"]:
        for rep in range(3):
            label = block.get("reference", "3_ns")
            source = next(r for r in reference["records"] + reference["late_block_records"]
                          if r["checkpoint"] == label and r["replicate"] == rep)
            record = {field: source[field] for field in blocks.REPRODUCTION_FIELDS}
            record.update(block=block["label"], replicate=rep,
                          n_raw=(block["stop"] - block["start"]) * 20,
                          n_retained=source["n_independent"])
            records.append(record)
    return records


def test_reproduction_checks_metrics_and_retained_counts():
    _, _, design = blocks.load_design()
    reference = _reference()
    records = _records_for_design(design, reference)
    blocks.check_reproduction(records, reference, design)
    records[0]["dg_kcal"] += .1
    with pytest.raises(ValueError, match="reference reproduction failed"):
        blocks.check_reproduction(records, reference, design)
    records[0]["dg_kcal"] -= .1
    records[0]["n_retained"] -= 1
    with pytest.raises(ValueError, match="sample counts"):
        blocks.check_reproduction(records, reference, design)


def test_paired_contrasts_and_repeat_sem_are_descriptive():
    _, _, design = blocks.load_design()
    records = _records_for_design(design, _reference())
    summaries, contrasts, contrast_summaries = blocks.summaries_and_contrasts(records, design)
    assert len(summaries) == 5
    assert len(contrasts) == 9
    assert len(contrast_summaries) == 3
    assert all(r["raw_columns_disjoint"] for r in contrasts)
    assert summaries[0]["mean_dg_kcal"] == pytest.approx(1.7042390759607582)
    assert summaries[2]["mean_dg_kcal"] == pytest.approx(.5870698195459618)
    assert summaries[2]["replicate_sem_kcal"] == pytest.approx(.10032918668852947)


def _mock_completed_pilot(tmp_path, monkeypatch):
    cfg, pilot, design = blocks.load_design()
    pilot = {**pilot, "_output_root": str(tmp_path)}
    reference = _reference()
    fep = tmp_path / "production/fep"
    folder = fep / "F64A/folded"
    folder.mkdir(parents=True)
    for window in range(20):
        for rep in range(3):
            (folder / f"w{window}_r{rep}.npz").write_bytes(b"mock validated input")
    analysis = tmp_path / "production/analysis"
    analysis.mkdir()
    (analysis / "f64a_extension.json").write_text(json.dumps(reference))
    (tmp_path / "production/stage_manifest.json").write_text(json.dumps(reference["code_identity"]))
    records = _records_for_design(design, reference)
    monkeypatch.setattr(blocks, "load_design", lambda config: (cfg, pilot, design))
    monkeypatch.setattr(blocks, "validate_production_outputs",
                        lambda config: (fep, {"code_identity": reference["code_identity"]}))
    monkeypatch.setattr(blocks, "_analysis_identity", lambda design: {"git_commit": "later_cpu"})
    monkeypatch.setattr(blocks, "solve_block", lambda cfg, pilot, fep, block, rep: copy.deepcopy(
        next(r for r in records if r["block"] == block["label"] and r["replicate"] == rep)))
    return fep, analysis


def test_cpu_pipeline_writes_new_report_preserves_primary_and_refuses_overwrite(tmp_path,
                                                                             monkeypatch):
    fep, analysis = _mock_completed_pilot(tmp_path, monkeypatch)
    primary = analysis / "f64a_extension.json"
    before = primary.read_bytes()
    output = blocks.analyze_blocks()
    result = json.loads(output.read_text())
    assert result["reproduction_passed"] is True
    assert len(result["records"]) == 15
    assert len(result["contrasts"]) == 9
    assert len([p for p in result["input_sha256"] if p.endswith(".npz")]) == 60
    assert result["simulation_code_identity"] == _reference()["code_identity"]
    assert result["analysis_identity"]["git_commit"] == "later_cpu"
    assert primary.read_bytes() == before
    with pytest.raises(FileExistsError, match="overwrite"):
        blocks.analyze_blocks()
    assert len(list(fep.rglob("*.npz"))) == 60


def test_validation_failure_precedes_solving_or_output(tmp_path, monkeypatch):
    _, analysis = _mock_completed_pilot(tmp_path, monkeypatch)

    def invalid(config):
        raise ValueError("mixed protocol")

    monkeypatch.setattr(blocks, "validate_production_outputs", invalid)
    monkeypatch.setattr(blocks, "solve_block", lambda *args: pytest.fail("must not solve"))
    with pytest.raises(ValueError, match="mixed protocol"):
        blocks.analyze_blocks()
    assert not (analysis / "f64a_blocks_v1.json").exists()


def test_primary_mismatch_is_fatal_without_partial_result(tmp_path, monkeypatch):
    _, analysis = _mock_completed_pilot(tmp_path, monkeypatch)
    (analysis / "f64a_extension.json").write_text("{}")
    with pytest.raises(ValueError, match="primary report differs"):
        blocks.analyze_blocks()
    assert not (analysis / "f64a_blocks_v1.json").exists()


def test_failed_fit_is_not_dropped_or_written_as_partial_result(tmp_path, monkeypatch):
    _, analysis = _mock_completed_pilot(tmp_path, monkeypatch)

    def failed(*args):
        raise FloatingPointError("failed synthetic MBAR fit")

    monkeypatch.setattr(blocks, "solve_block", failed)
    with pytest.raises(FloatingPointError, match="failed synthetic"):
        blocks.analyze_blocks()
    assert not (analysis / "f64a_blocks_v1.json").exists()


def test_input_changes_during_fits_are_fatal(tmp_path, monkeypatch):
    fep, analysis = _mock_completed_pilot(tmp_path, monkeypatch)
    original_solve = blocks.solve_block

    def changed(*args):
        result = original_solve(*args)
        (fep / "F64A/folded/w0_r0.npz").write_bytes(b"changed mock input")
        return result

    monkeypatch.setattr(blocks, "solve_block", changed)
    with pytest.raises(ValueError, match="input changed"):
        blocks.analyze_blocks()
    assert not (analysis / "f64a_blocks_v1.json").exists()


def test_frozen_report_and_gpu_configuration_still_match():
    reference = _reference()
    assert reference["code_identity"]["git_commit"] == "80097dcb2ff90ebfdedfadb6e4d9a866a3d9b55a"
    assert reference["pilot_config_sha256"] == sha256_file(
        ROOT / "pilots/f64a_sampling_extension_v1/config.yaml")
    assert reference["base_config_sha256"] == sha256_file(ROOT / "config/pipeline.yaml")
    assert reference["protocol"] == "bf6841ccb3b9de79"


def test_snakemake_cpu_target_uses_completed_inputs_and_never_schedules_gpu(tmp_path):
    executable = shutil.which("snakemake")
    if executable is None:
        pytest.skip("Snakemake not installed")
    (tmp_path / "config").mkdir()
    (tmp_path / "config/pipeline.yaml").write_bytes((ROOT / "config/pipeline.yaml").read_bytes())
    (tmp_path / "workflow").mkdir()
    (tmp_path / "workflow/Snakefile").write_bytes((ROOT / "workflow/Snakefile").read_bytes())
    package = tmp_path / "pilots/f64a_sampling_extension_v1"
    package.mkdir(parents=True)
    pilot = yaml.safe_load((ROOT / "pilots/f64a_sampling_extension_v1/config.yaml").read_text())
    pilot["output_root"] = "mock_completed"
    (package / "config.yaml").write_text(yaml.safe_dump(pilot))
    (package / "cpu_blocks.yaml").write_bytes(blocks.DEFAULT_CONFIG.read_bytes())
    (tmp_path / "src").symlink_to(ROOT / "src", target_is_directory=True)
    (tmp_path / "docs").symlink_to(ROOT / "docs", target_is_directory=True)
    (tmp_path / "data").symlink_to(ROOT / "data", target_is_directory=True)
    completed = tmp_path / "mock_completed/production"
    (completed / "analysis").mkdir(parents=True)
    (completed / "stage_manifest.json").write_text("{}")
    (completed / "analysis/f64a_extension.json").write_text("{}")
    folder = completed / "fep/F64A/folded"
    folder.mkdir(parents=True)
    for window in range(20):
        for rep in range(3):
            (folder / f"w{window}_r{rep}.npz").write_bytes(b"mock completed input")
    command = [executable, "f64a_extension_blocks", "-s", "workflow/Snakefile", "-n",
               "--cores", "1", "--runtime-source-cache-path", str(tmp_path / "cache"),
               "--allowed-rules", "f64a_extension_blocks"]
    dry_run = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    output = dry_run.stdout + dry_run.stderr
    assert dry_run.returncode == 0, output
    assert "rule f64a_extension_blocks:" in output
    assert "rule f64a_extension_window:" not in output
    assert "rule f64a_extension_analysis:" not in output
    (folder / "w0_r0.npz").unlink()
    missing = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    output = missing.stdout + missing.stderr
    assert missing.returncode != 0
    assert "will not schedule GPU work" in output
    assert "rule f64a_extension_window:" not in output
