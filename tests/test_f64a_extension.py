"""Scientific and isolation contracts for the F64A exact-path pilot."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from src.fep.f64a_extension import (
    _extension_protocol,
    _validate_output_npz,
    _validate_baseline_mdp,
    check_first_light,
    analyze,
    load_pilot,
    reconcile_energy_history,
    run_task,
    sha256_file,
    stage,
)

ROOT = Path(__file__).resolve().parents[1]


def _identity() -> dict:
    return {
        "git_commit": "a" * 40,
        "pilot_code_sha256": "b" * 64,
        "pilot_paths_clean": True,
        "pilot_git_status": [],
    }


def _pilot_config(tmp_path: Path) -> Path:
    cfg = yaml.safe_load(
        (ROOT / "pilots/f64a_sampling_extension_v1/config.yaml").read_text()
    )
    cfg["base_config"] = str(ROOT / "config/pipeline.yaml")
    cfg["output_root"] = str(tmp_path / "isolated-output")
    path = tmp_path / "pilot.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def _source_task(root: Path, window: int = 0, rep: int = 0) -> None:
    name = f"w{window}_r{rep}"
    run = root / name
    run.mkdir(parents=True)
    np.savez(
        root / f"{name}.npz",
        lambda_index=window,
        n_states=20,
        u_kn_window=np.zeros((20, 3001)),
        provenance="gromacs_pmx",
        protocol="822108e9db71124d",
    )
    (run / "prod.mdp").write_text(
        "dt = 0.002\n"
        "nsteps = 1750000\n"
        "nstdhdl = 500\n"
        "nstxout-compressed = 0\n"
        f"init-lambda-state = {window}\n"
    )
    (run / "protocol.sha").write_text("822108e9db71124d\n")
    for filename in ("prod.tpr", "prod.cpt", "prod.log", "prod.edr", "dhdl.xvg"):
        (run / filename).write_bytes(f"test-{filename}".encode())


def test_registered_pilot_is_isolated_from_historical_results():
    _, pilot = load_pilot()
    output = Path(pilot["_output_root"])
    historical = (ROOT / "results/fep").resolve()
    assert output != historical
    assert historical not in output.parents
    assert pilot["variant"] == "F64A"
    assert pilot["leg"] == "folded"
    assert pilot["lambda_windows"] * pilot["replicates"] == 60


def test_load_pilot_rejects_scientific_override(tmp_path):
    path = _pilot_config(tmp_path)
    cfg = yaml.safe_load(path.read_text())
    cfg["extension_ns"] = 5.0
    path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError, match="registered contract"):
        load_pilot(path)


def test_stage_first_light_copies_and_hashes_without_touching_source(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    source = tmp_path / "archive" / "folded"
    _source_task(source)
    before = {path.relative_to(source): path.read_bytes() for path in source.rglob("*")
              if path.is_file()}
    config = _pilot_config(tmp_path)
    monkeypatch.setattr(extension, "validate_source",
                        lambda *args, **kwargs: {"tasks": 1, "files": 8, "bytes": 1})
    monkeypatch.setattr(extension, "_require_committed_code", _identity)

    root = stage(source, config, first_light=True)

    assert (root / "stage_manifest.json").is_file()
    assert (root / "source/F64A/folded/w0_r0.npz").is_file()
    after = {path.relative_to(source): path.read_bytes() for path in source.rglob("*")
             if path.is_file()}
    assert after == before
    # Idempotent re-entry verifies the completed stage rather than overwriting it.
    assert stage(source, config, first_light=True) == root


def test_stage_preflights_missing_checkpoint_without_creating_destination(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    source = tmp_path / "archive" / "folded"
    _source_task(source)
    (source / "w0_r0/prod.cpt").unlink()
    config = _pilot_config(tmp_path)
    monkeypatch.setattr(extension, "_require_committed_code", _identity)
    with pytest.raises(FileNotFoundError, match="source is incomplete"):
        stage(source, config, first_light=True)
    assert not (tmp_path / "isolated-output/first_light").exists()


def test_baseline_mdp_contract_rejects_wrong_lambda(tmp_path):
    path = tmp_path / "prod.mdp"
    path.write_text(
        "dt = 0.002\nnsteps = 1750000\nnstdhdl = 500\n"
        "nstxout-compressed = 0\ninit-lambda-state = 1\n"
    )
    with pytest.raises(ValueError, match="continuation contract"):
        _validate_baseline_mdp(path, window=0)


def test_extension_protocol_distinguishes_first_light_from_production():
    _, pilot = load_pilot()
    assert _extension_protocol(pilot, 10.0) != _extension_protocol(pilot, 6000.0)


def test_run_task_appends_in_private_copy_and_combines_3_plus_6_ns(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    config = _pilot_config(tmp_path)
    _, pilot = load_pilot(config)
    root = Path(pilot["_output_root"]) / "production"
    source = root / "source/F64A/folded"
    _source_task(source)
    identity = _identity()
    (root / "stage_manifest.json").parent.mkdir(parents=True, exist_ok=True)
    (root / "stage_manifest.json").write_text(json.dumps({
        "pilot_id": pilot["pilot_id"],
        "pilot_config_sha256": sha256_file(config),
        "base_config_sha256": sha256_file(pilot["_base_config_path"]),
        "code_identity": identity,
        "files": [{
            "path": "source/F64A/folded/w0_r0.npz",
            "sha256": sha256_file(source / "w0_r0.npz"),
        }],
    }))

    def fake_run(argv, cwd, **kwargs):
        if "convert-tpr" in argv:
            requested = Path(argv[argv.index("-o") + 1])
            assert requested.is_absolute()
            assert requested.suffix == ".tpr"
            requested.write_bytes(b"extended tpr")
        return ""

    monkeypatch.setattr(extension, "_require_committed_code", _identity)
    monkeypatch.setattr(extension, "gmx_command", lambda cfg: "gmx")
    monkeypatch.setattr(extension, "_run", fake_run)
    monkeypatch.setattr(extension, "_run_extension_mdrun", lambda cfg, run_dir: None)
    monkeypatch.setattr(extension, "reconcile_energy_history",
                        lambda *args, **kwargs: np.ones((20, 6000)))

    out = run_task(0, 0, config)

    with np.load(out) as data:
        assert data["u_kn_window"].shape == (20, 9001)
        assert int(data["baseline_samples"]) == 3001
        assert int(data["extension_samples"]) == 6000
        assert str(data["git_commit"]) == identity["git_commit"]
        assert str(data["pilot_code_sha256"]) == identity["pilot_code_sha256"]
    private_run = root / "fep/F64A/folded/w0_r0"
    assert (private_run / "prod.cpt").read_bytes() == b"test-prod.cpt"
    assert (private_run / "prod.log").read_bytes() == b"test-prod.log"
    assert (private_run / "prod.edr").read_bytes() == b"test-prod.edr"


def test_submission_geometry_and_gpu_pin_are_fixed():
    text = (ROOT / "pilots/f64a_sampling_extension_v1/submit_array.sh").read_text()
    assert "#$ -t 1-60" in text
    assert "#$ -tc 8" in text
    assert "#$ -l gpu_type=L40S" in text
    assert "#$ -l h_rt=12:00:00" in text
    assert "check-first-light" in text
    assert "results/fep/F64A" not in text
    first_light = (ROOT / "pilots/f64a_sampling_extension_v1/submit_first_light.sh").read_text()
    assert "SGE_TASK_ID" not in first_light
    snakefile = (ROOT / "workflow/Snakefile").read_text()
    assert "rule f64a_extension_window:" in snakefile
    assert "rule f64a_extension_analysis:" in snakefile


def test_first_light_gate_requires_finite_expected_output(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    config = _pilot_config(tmp_path)
    _, pilot = load_pilot(config)
    identity = _identity()
    root = Path(pilot["_output_root"]) / "first_light"
    path = (Path(pilot["_output_root"]) / "first_light/fep/F64A/folded/w0_r0.npz")
    path.parent.mkdir(parents=True)
    source_hash = "c" * 64
    manifest = {
        "pilot_id": pilot["pilot_id"],
        "pilot_config_sha256": sha256_file(config),
        "base_config_sha256": sha256_file(pilot["_base_config_path"]),
        "code_identity": identity,
        "files": [{"path": "source/F64A/folded/w0_r0.npz", "sha256": source_hash}],
    }
    (root / "stage_manifest.json").write_text(json.dumps(manifest))
    np.savez(
        path,
        lambda_index=0,
        n_states=20,
        u_kn_window=np.zeros((20, 3011)),
        provenance="gromacs_pmx",
        pilot_id=pilot["pilot_id"],
        protocol=_extension_protocol(pilot, 10.0),
        source_protocol=pilot["source_protocol"],
        source_npz_sha256=source_hash,
        baseline_samples=3001,
        extension_samples=10,
        git_commit=identity["git_commit"],
        pilot_code_sha256=identity["pilot_code_sha256"],
        pilot_config_sha256=sha256_file(config),
        base_config_sha256=sha256_file(pilot["_base_config_path"]),
    )
    monkeypatch.setattr(extension, "_require_committed_code", _identity)
    assert check_first_light(config) == path
    path.unlink()
    with pytest.raises(FileNotFoundError, match="first-light result is missing"):
        check_first_light(config)


def test_reconcile_energy_history_requires_exact_grid_and_baseline(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    config = _pilot_config(tmp_path)
    cfg, pilot = load_pilot(config)
    baseline = np.arange(20 * 3001, dtype=float).reshape(20, 3001)
    npz = tmp_path / "w0_r0.npz"
    np.savez(npz, u_kn_window=baseline)
    xvg = tmp_path / "dhdl.xvg"
    xvg.write_text("\n".join(f"{i}.0000 0" for i in range(3501)) + "\n")
    full = np.zeros((20, 3501))
    full[:, 500:] = baseline
    monkeypatch.setattr(extension, "dhdl_to_u_kn", lambda *args, **kwargs: full)
    assert reconcile_energy_history(xvg, npz, cfg, pilot, 3500.0).shape == (20, 0)

    lines = xvg.read_text().splitlines()
    lines[2000] = "2001.0000 0"
    xvg.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="time grid is not exact"):
        reconcile_energy_history(xvg, npz, cfg, pilot, 3500.0)


def test_existing_output_rejects_wrong_protocol(tmp_path):
    config = _pilot_config(tmp_path)
    _, pilot = load_pilot(config)
    path = tmp_path / "w0_r0.npz"
    identity = _identity()
    manifest = {
        "pilot_config_sha256": sha256_file(config),
        "base_config_sha256": sha256_file(pilot["_base_config_path"]),
        "code_identity": identity,
        "files": [{"path": "source/F64A/folded/w0_r0.npz", "sha256": "c" * 64}],
    }
    np.savez(
        path,
        lambda_index=0, n_states=20, u_kn_window=np.zeros((20, 9001)),
        provenance="gromacs_pmx", protocol="wrong", pilot_id=pilot["pilot_id"],
        source_protocol=pilot["source_protocol"], source_npz_sha256="c" * 64,
        baseline_samples=3001, extension_samples=6000,
        git_commit=identity["git_commit"], pilot_code_sha256=identity["pilot_code_sha256"],
        pilot_config_sha256=sha256_file(config),
        base_config_sha256=sha256_file(pilot["_base_config_path"]),
    )
    with pytest.raises(ValueError, match="protocol=wrong"):
        _validate_output_npz(path, pilot, 6000.0, manifest)


def test_analyze_runs_complete_input_validation_before_mbar(tmp_path, monkeypatch):
    import src.fep.f64a_extension as extension

    config = _pilot_config(tmp_path)
    monkeypatch.setattr(
        extension, "validate_production_outputs",
        lambda config: (_ for _ in ()).throw(ValueError("mixed protocol test")),
    )
    with pytest.raises(ValueError, match="mixed protocol test"):
        analyze(config)
