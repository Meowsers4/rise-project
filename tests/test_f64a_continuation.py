"""Synthetic tests of v2 isolation, admission, boundary and CPU-analysis contracts."""

from __future__ import annotations

import copy
import json
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from src.analysis import f64a_continuation as analysis
from src.fep import f64a_continuation as continuation
from src.fep.f64a_extension import ROOT, sha256_file


def _identity():
    return {"git_commit": "synthetic", "source_sha256": {"test": "synthetic"}}


def _config(tmp_path, monkeypatch):
    values = yaml.safe_load(continuation.DEFAULT_CONFIG.read_text())
    old = yaml.safe_load((ROOT / values["source_config"]).read_text())
    old["base_config"] = str(ROOT / "config/pipeline.yaml")
    old["output_root"] = str(tmp_path / "results/pilots/v1")
    old_path = tmp_path / "v1_config.yaml"
    old_path.write_text(yaml.safe_dump(old))
    values.update(base_config=str(ROOT / "config/pipeline.yaml"), source_config=str(old_path),
                  source_production=str(tmp_path / "results/pilots/v1/production"),
                  output_root=str(tmp_path / "results/pilots/v2"))
    config = tmp_path / "v2_config.yaml"
    config.write_text(yaml.safe_dump(values))
    monkeypatch.setattr(continuation, "ROOT", tmp_path)
    monkeypatch.setattr(continuation, "code_identity", lambda pilot: _identity())
    return config


def _accounting(task_count=60, wall=1000):
    return "\n".join(f"{'=' * 30}\njobnumber 7589742\ntaskid {i}\nfailed 0\n"
                     f"exit_status 0\nru_wallclock {wall}\ncpu {wall * 8}\n"
                     for i in range(1, task_count + 1))


def test_selected_config_fixed_science_and_distinct_roots():
    _, pilot = continuation.load_pilot()
    assert pilot["target_samples_per_window"] == 15001
    assert pilot["source_protocol"] == "bf6841ccb3b9de79"
    assert pilot["_source_root"] not in pilot["_output_root"]


@pytest.mark.parametrize("change", [{"target_ns": 18.}, {"replicates": 5},
                                    {"gpu_type": "Blackwell"}, {"output_root": "results/fep"}])
def test_invalid_contract_or_output_refused(tmp_path, monkeypatch, change):
    config = _config(tmp_path, monkeypatch)
    values = yaml.safe_load(config.read_text())
    values.update(change)
    config.write_text(yaml.safe_dump(values))
    with pytest.raises(ValueError):
        continuation.load_pilot(config)


@pytest.mark.parametrize("corrupt", ["missing", "failed", "wrong_job", "duplicate", "cpu_nan"])
def test_accounting_refuses_incomplete_or_failed_inputs(corrupt):
    value = _accounting()
    if corrupt == "missing":
        value = _accounting(59)
    elif corrupt == "failed":
        value = value.replace("failed 0", "failed 1", 1)
    elif corrupt == "wrong_job":
        value = value.replace("7589742", "123", 1)
    elif corrupt == "duplicate":
        value = value.replace("taskid 60", "taskid 59")
    else:
        value = value.replace("cpu 8000", "cpu nan", 1)
    with pytest.raises(ValueError):
        continuation.parse_accounting(value, 7589742)


def test_accounting_uses_occupied_runtime_not_queue_elapsed():
    value = continuation.parse_accounting(_accounting(), 7589742)
    assert value["gpu_hours"] == pytest.approx(60 * 1000 / 3600)
    assert value["maximum_task_seconds"] == 1000


def _source_history(tmp_path, monkeypatch, end=9510.):
    cfg, pilot = continuation.load_pilot()
    npz = tmp_path / "source.npz"
    frozen = np.zeros((20, 9001))
    np.savez(npz, u_kn_window=frozen, append_boundary_max_abs_reduced_potential_difference=.2)
    times = np.arange(0., end + 1)
    parsed = np.zeros((20, len(times)))
    parsed[0, 3500] = .2
    if end > 9500:
        parsed[0, 9500] = .3
    xvg = tmp_path / "dhdl.xvg"
    monkeypatch.setattr(continuation.v1, "_xvg_times", lambda path: times)
    monkeypatch.setattr(continuation, "dhdl_to_u_kn", lambda *args: parsed)
    return cfg, pilot, npz, xvg, parsed, times


def test_two_boundaries_audited_and_never_analyzed_as_new_columns(tmp_path, monkeypatch):
    cfg, pilot, npz, xvg, parsed, _ = _source_history(tmp_path, monkeypatch)
    parsed[:, 9501:] = 4.
    new, deltas = continuation.reconcile(xvg, npz, cfg, pilot, 9510.)
    assert new.shape == (20, 10)
    assert np.all(new == 4.)
    assert deltas == {"3500": .2, "9500": .3}


@pytest.mark.parametrize("corrupt", ["old_delta", "history", "duplicate", "nonfinite", "source_endpoint"])
def test_energy_integrity_refuses_every_unregistered_change(tmp_path, monkeypatch, corrupt):
    end = 9500. if corrupt == "source_endpoint" else 9510.
    cfg, pilot, npz, xvg, parsed, times = _source_history(tmp_path, monkeypatch, end)
    if corrupt == "old_delta":
        parsed[1, 3500] = .5
    elif corrupt == "history":
        parsed[0, 9499] = .01
    elif corrupt == "duplicate":
        times[-1] = times[-2]
    elif corrupt == "nonfinite":
        parsed[0, -1] = np.nan
    else:
        parsed[0, 9500] = .01
    with pytest.raises(ValueError):
        continuation.reconcile(xvg, npz, cfg, pilot, end)


def test_real_tpr_fields_and_checkpoint_time_are_checked(tmp_path, monkeypatch):
    cfg, _ = continuation.load_pilot()
    dump = ("inputrec:\ndt = 0.002\nnsteps = 4750000\ninit-step = 0\ntinit = 0\n"
            "init-lambda-state = 7\nnstdhdl = 500\nnstxout = 0\nnstvout = 0\n"
            "nstfout = 0\nnstxout-compressed = 0\n")
    monkeypatch.setattr(continuation, "gmx_command", lambda cfg: "gmx")
    monkeypatch.setattr(continuation, "dump_tpr", lambda *a, **k: dump)
    continuation.check_tpr(cfg, tmp_path / "extension.tpr", 7, 9500.)
    with pytest.raises(ValueError):
        continuation.check_tpr(cfg, tmp_path / "prod.tpr", 7, 15500.)
    monkeypatch.setattr(continuation, "_run", lambda *a, **k: "t = 9500\n")
    continuation.checkpoint_time(cfg, tmp_path / "prod.cpt", 9500.)
    with pytest.raises(ValueError):
        continuation.checkpoint_time(cfg, tmp_path / "prod.cpt", 3500.)


def test_tpr_model_fingerprint_allows_only_duration_change():
    original = "source.tpr\ninputrec:\n   nsteps = 4750000\n   dt = 0.002\ntopology:\nx = 0.1\n"
    extended = original.replace("source.tpr", "extension.tpr").replace("4750000", "7750000")
    assert continuation.tpr_model_hash(original) == continuation.tpr_model_hash(extended)
    for modified in (extended.replace("0.002", "0.004"), extended.replace("x = 0.1", "x = 0.2")):
        assert continuation.tpr_model_hash(original) != continuation.tpr_model_hash(modified)


def _admission_fixture(tmp_path, monkeypatch):
    config = _config(tmp_path, monkeypatch)
    files = []
    source = tmp_path / "source/F64A/folded"
    run = source / "w0_r0"
    run.mkdir(parents=True)
    np.savez(source / "w0_r0.npz", u_kn_window=np.zeros((20, 9001)),
             append_boundary_max_abs_reduced_potential_difference=0.)
    for name in ("extension.tpr", "prod.cpt", "prod.log", "prod.edr", "dhdl.xvg", "extension_protocol.sha"):
        (run / name).write_bytes(b"synthetic source file")
    for path, relative in continuation.source_files(source, 0, 0):
        files.append(dict(source=str(path), relative=str(relative), sha256=sha256_file(path),
                          size_bytes=path.stat().st_size))
    source_inventory = dict(files=files, tasks=60, bytes=sum(e["size_bytes"] for e in files),
                            source_code_identity={"git_commit": "v1_synthetic"},
                            tpr_checks=[dict(task='w0_r0', values={'model_sha256_excluding_nsteps': 'synthetic'})])
    monkeypatch.setattr(continuation, "validate_source", lambda config: source_inventory)
    evidence = []
    for name, content in (("accounting", _accounting()), ("quota", "synthetic quota"),
                          ("balance", "synthetic balance")):
        path = tmp_path / (name + ".txt")
        path.write_text(content)
        evidence.append(path)
    return config, evidence, source


def test_admission_budget_is_a_finite_reservation_and_no_overwrite(tmp_path, monkeypatch):
    config, evidence, _ = _admission_fixture(tmp_path, monkeypatch)
    path = continuation.admit(config, *evidence, available_gib=10., su_balance=10000.,
                             su_rate=1., gpu_hour_cap=100., su_cap=800., max_attempts=2)
    value = json.loads(path.read_text())
    assert value["maximum_reserved_gpu_hours"] <= 100.
    assert value["maximum_reserved_su"] <= 800.
    assert value["task_seconds"] <= 43200
    with pytest.raises(FileExistsError):
        continuation.admit(config, *evidence, 10., 10000., 1., 100., 800., 2)


@pytest.mark.parametrize("values", [(10., 10000., 1., 5., 800., 2),
                                    (.00001, 10000., 1., 100., 800., 2),
                                    (10., 10000., 1., 100., 10., 2)])
def test_admission_rejects_runtime_quota_or_su_without_approval(tmp_path, monkeypatch, values):
    config, evidence, _ = _admission_fixture(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        continuation.admit(config, *evidence, *values)
    _, pilot = continuation.load_pilot(config)
    assert not (Path(pilot["_output_root"]) / "admission.json").exists()


def _staged_light(tmp_path, monkeypatch):
    config, evidence, source = _admission_fixture(tmp_path, monkeypatch)
    continuation.admit(config, *evidence, 10., 10000., 1., 100., 800., 2)
    root = continuation.stage(config, first_light=True)
    return config, root, source


def test_stage_copies_all_sources_without_mutation_and_guards_first_light(tmp_path, monkeypatch):
    config, root, source = _staged_light(tmp_path, monkeypatch)
    hashes = {p: sha256_file(p) for p in source.rglob('*') if p.is_file()}
    manifest = json.loads((root / "stage_manifest.json").read_text())
    assert len(manifest["files"]) == 7
    assert continuation.stage(config, True) == root
    assert hashes == {p: sha256_file(p) for p in hashes}
    with pytest.raises(FileNotFoundError):
        continuation.stage(config, False)
    (root / "source/F64A/folded/w0_r0/prod.cpt").write_bytes(b"corrupt")
    with pytest.raises(ValueError):
        continuation.stage(config, True)


def _run_light(tmp_path, monkeypatch, returncode=0, output=""):
    config, root, source = _staged_light(tmp_path, monkeypatch)
    monkeypatch.setenv("JOB_ID", "123")
    (root.parent / "first_light_submission.json").write_text(json.dumps({"job_id": "123"}))
    monkeypatch.setattr(continuation, "gmx_command", lambda cfg: "gmx")
    seen = []

    def convert(argv, **kwargs):
        assert "convert-tpr" in argv and argv[argv.index('-s') + 1].endswith('/extension.tpr')
        Path(argv[argv.index('-o') + 1]).write_bytes(b"synthetic extended TPR")
        seen.append(argv)

    monkeypatch.setattr(continuation, "_run", convert)
    monkeypatch.setattr(continuation, "check_tpr", lambda *a: {})
    monkeypatch.setattr(continuation, "checkpoint_time", lambda *a: None)
    monkeypatch.setattr(continuation, "mdrun_argv", lambda cfg: ["gmx", "mdrun"])
    monkeypatch.setattr(continuation.subprocess, "run", lambda argv, **k: SimpleNamespace(
        returncode=returncode, stdout=output, stderr=""))
    monkeypatch.setattr(continuation, "reconcile", lambda *a: (np.ones((20, 10)), {"3500": 0., "9500": .3}))
    return config, root, source, seen


def test_private_run_uses_9ns_tpr_preserves_prefix_and_validates_existing_result(tmp_path, monkeypatch):
    config, root, source, seen = _run_light(tmp_path, monkeypatch)
    before = {p: sha256_file(p) for p in source.rglob('*') if p.is_file()}
    output = continuation.run_task(0, 0, config, first_light=True)
    assert seen
    with np.load(output) as data:
        assert data['u_kn_window'].shape == (20, 9011)
        assert np.all(data['u_kn_window'][:, :9001] == 0.)
        assert np.all(data['u_kn_window'][:, 9001:] == 1.)
    assert continuation.run_task(0, 0, config, True) == output
    assert continuation.check_first_light(config) == output
    assert before == {p: sha256_file(p) for p in before}
    with np.load(output) as data:
        changed = {k: data[k] for k in data.files}
    changed['protocol'] = 'wrong'
    np.savez(output, **changed)
    with pytest.raises(ValueError, match="protocol"):
        continuation.run_task(0, 0, config, True)


def test_only_gpu_placement_failure_reschedules_and_attempt_budget_is_finite(tmp_path, monkeypatch):
    config, root, _, _ = _run_light(tmp_path, monkeypatch, 1, "no GPU is detected")
    with pytest.raises(continuation.GpuUnavailableError):
        continuation.run_task(0, 0, config, True)
    with pytest.raises(continuation.ToolError):
        continuation.run_task(0, 0, config, True)
    with pytest.raises(continuation.ToolError, match="exhausted"):
        continuation.run_task(0, 0, config, True)
    assert len(json.loads((root / 'fep/F64A/folded/w0_r0/attempts.json').read_text())) == 2


def test_nonplacement_failure_does_not_reschedule(tmp_path, monkeypatch):
    config, _, _, _ = _run_light(tmp_path, monkeypatch, 1, 'LINCS warning cudaErrorInvalidPtx')
    with pytest.raises(continuation.ToolError):
        continuation.run_task(0, 0, config, True)


def test_analysis_design_has_54_fits_and_disjoint_contrasts():
    _, pilot = continuation.load_pilot()
    design = yaml.safe_load(continuation.repo_path(pilot['analysis_config']).read_text())
    assert 3 * 3 * len(design['blocks']) + 3 * len(design['cumulative']) == 54
    by_label = {b['label']: b for b in design['blocks']}
    for c in design['contrasts']:
        assert by_label[c['earlier']]['stop'] <= by_label[c['later']]['start']
    assert by_label['block_13p5_15_ns']['stop'] == 15001


@pytest.mark.parametrize('corruption', ['path', 'policy', 'slice', 'overlap', 'count'])
def test_analysis_design_refuses_unsafe_or_incomplete_plan(tmp_path, corruption):
    _, pilot = continuation.load_pilot()
    path = continuation.repo_path(pilot['analysis_config'])
    design = analysis.load_design(path, pilot)
    if corruption == 'path':
        design['output_name'] = '../outside.json'
    elif corruption == 'policy':
        design['policies'][2]['trim_fraction'] = .5
    elif corruption == 'slice':
        design['blocks'][0]['stop'] = 20000
    elif corruption == 'overlap':
        design['blocks'][2]['start'] = 9002
    else:
        design['cumulative'].pop()
    path = tmp_path / 'analysis.yaml'
    path.write_text(yaml.safe_dump(design))
    with pytest.raises(ValueError):
        analysis.load_design(path, pilot)


def test_v1_reproduction_includes_all_six_metrics_and_counts():
    reference = json.loads((ROOT / 'docs/f64a_extension_review/f64a_extension.json').read_text())
    records = []
    for label, checkpoint in [('cumulative_9_ns', '9_ns'), ('block_6_9_ns', 'final_3_ns')]:
        for r in reference['records'] + reference['late_block_records']:
            if r['checkpoint'] == checkpoint:
                records.append({**r, 'policy': 'adaptive', 'block': label, 'n_retained': r['n_independent']})
    analysis.check_reproduction(records, reference, {'reproduction_atol': 1e-8})
    records[0]['dg_kcal'] += .1
    with pytest.raises(ValueError):
        analysis.check_reproduction(records, reference, {'reproduction_atol': 1e-8})


def test_existing_v1_configs_and_source_are_unchanged():
    report = json.loads((ROOT / 'docs/f64a_extension_review/f64a_selection_sensitivity_v1.json').read_text())
    for path, digest in report['analysis_identity']['source_sha256'].items():
        assert sha256_file(ROOT / path) == digest


@pytest.mark.parametrize('completed', [False, True])
def test_cpu_analysis_only_consumes_completed_inputs(tmp_path, completed):
    executable = shutil.which('snakemake')
    if not executable:
        pytest.skip('Snakemake unavailable')
    for name in ('config', 'workflow', 'pilots/f64a_sampling_extension_v1', 'pilots/f64a_sampling_extension_v2'):
        (tmp_path / name).mkdir(parents=True)
    for name in ('config/pipeline.yaml', 'workflow/Snakefile',
                 'pilots/f64a_sampling_extension_v1/config.yaml',
                 'pilots/f64a_sampling_extension_v1/cpu_blocks.yaml',
                 'pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml',
                 'pilots/f64a_sampling_extension_v2/config.yaml',
                 'pilots/f64a_sampling_extension_v2/PREREGISTRATION.md',
                 'pilots/f64a_sampling_extension_v2/analysis.yaml'):
        (tmp_path / name).write_bytes((ROOT / name).read_bytes())
    for name in ('src', 'docs', 'data'):
        (tmp_path / name).symlink_to(ROOT / name, target_is_directory=True)
    if completed:
        root = tmp_path / 'results/pilots/f64a_sampling_extension_v2'
        folder = root / 'production/fep/F64A/folded'
        folder.mkdir(parents=True)
        (root / 'admission.json').write_text('{}')
        (root / 'production/stage_manifest.json').write_text('{}')
        for w in range(20):
            for r in range(3):
                (folder / f'w{w}_r{r}.npz').write_bytes(b'synthetic completed input')
    command = [executable, 'f64a_v2_analysis', '-s', 'workflow/Snakefile', '-n', '--cores', '1',
               '--runtime-source-cache-path', str(tmp_path / 'cache'), '--allowed-rules', 'f64a_v2_analysis']
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    text = result.stdout + result.stderr
    if completed:
        assert result.returncode == 0, text
        assert 'rule f64a_v2_analysis:' in text
    else:
        assert result.returncode != 0
        assert 'will not schedule GPU work' in text
    assert 'rule f64a_v2_window:' not in text
    assert 'rule f64a_v2_first_light:' not in text


def _mock_analysis(tmp_path, monkeypatch):
    cfg, pilot = continuation.load_pilot()
    pilot = {**pilot, '_output_root': str(tmp_path)}
    reference = json.loads((ROOT / pilot['reference_primary']).read_text())
    folder = tmp_path / 'production/fep/F64A/folded'
    folder.mkdir(parents=True)
    for w in range(20):
        for r in range(3):
            (folder / f'w{w}_r{r}.npz').write_bytes(b'synthetic completed input')
    (tmp_path / 'production/stage_manifest.json').write_text('{}')
    (tmp_path / 'admission.json').write_text('{}')
    stage = dict(files=[], code_identity=_identity(), source_code_identity={'git_commit': 'v1'})
    monkeypatch.setattr(continuation, 'load_pilot', lambda config: (cfg, pilot))
    monkeypatch.setattr(continuation, 'validate_outputs', lambda config: (folder.parents[1], stage))
    monkeypatch.setattr(analysis.blocks, '_analysis_identity', lambda design: _identity())

    def solve(cfg, pilot, fep, block, rep, **kwargs):
        checkpoint = 'final_3_ns' if block['label'] == 'block_6_9_ns' else '9_ns'
        record = copy.deepcopy(next(r for r in reference['records'] + reference['late_block_records']
                                    if r['checkpoint'] == checkpoint and r['replicate'] == rep))
        record.update(block=block['label'], n_raw=(block['stop'] - block['start']) * 20,
                      n_retained=record['n_independent'],
                      adjacent_diagnostics=[dict(lower_window=k, upper_window=k+1,
                          signed_discrepancy_kcal=0., forward_dg_kcal=0., reverse_dg_kcal=0.)
                          for k in range(19)])
        return record

    monkeypatch.setattr(analysis.blocks, 'solve_block', solve)
    return folder


def test_cpu_pipeline_writes_all_54_fits_without_touching_inputs(tmp_path, monkeypatch):
    folder = _mock_analysis(tmp_path, monkeypatch)
    before = {p: p.read_bytes() for p in folder.iterdir()}
    result = json.loads(analysis.analyze().read_text())
    assert len(result['records']) == 54
    assert len(result['contrast_summaries']) == 9
    assert len(result['contrasts']) == 27
    assert len(result['local_changes']) == 513
    assert result['v1_reproduction_passed']
    assert before == {p: p.read_bytes() for p in before}
    with pytest.raises(FileExistsError):
        analysis.analyze()


@pytest.mark.parametrize('corruption', ['fit', 'reproduction', 'mutation'])
def test_cpu_pipeline_failure_never_writes_partial_result(tmp_path, monkeypatch, corruption):
    folder = _mock_analysis(tmp_path, monkeypatch)
    original = analysis.blocks.solve_block

    def solve(*args, **kwargs):
        if corruption == 'fit':
            raise ValueError('synthetic failed fit')
        record = original(*args, **kwargs)
        if corruption == 'reproduction':
            record['dg_kcal'] += .1
        else:
            (folder / 'w0_r0.npz').write_bytes(b'changed input')
        return record

    monkeypatch.setattr(analysis.blocks, 'solve_block', solve)
    with pytest.raises(ValueError):
        analysis.analyze()
    assert not (tmp_path / 'production/analysis/f64a_continuation_v2.json').exists()


def test_guarded_submission_caps_walltime_records_identity_and_refuses_duplicate(tmp_path, monkeypatch):
    config, root, _ = _staged_light(tmp_path, monkeypatch)
    seen = []

    def scheduler(argv, **kwargs):
        seen.append(argv)
        if argv[0] == 'qrls':
            assert (root.parent / 'first_light_submission.json').exists()
        return SimpleNamespace(stdout='456\n', stderr='', returncode=0)

    monkeypatch.setattr(continuation.subprocess, 'run', scheduler)
    assert continuation.submit(config, True) == '456'
    assert 'h_rt=01:00:00' in seen[0]
    assert '-h' in seen[0]
    assert seen[1] == ['qrls', '456']
    assert json.loads((root.parent / 'first_light_submission.json').read_text())['job_id'] == '456'
    with pytest.raises(FileExistsError):
        continuation.submit(config, True)
    assert len(seen) == 2


def test_retry_rejects_active_duplicate_and_exhausted_ledger(tmp_path, monkeypatch):
    config, root, _ = _staged_light(tmp_path, monkeypatch)
    (root.parent / 'first_light_submission.json').write_text(json.dumps({'job_id': '456'}))
    monkeypatch.setattr(continuation.subprocess, 'run', lambda argv, **kwargs: SimpleNamespace(
        stdout='456 0.1 running_job\n', stderr='', returncode=0))
    with pytest.raises(RuntimeError, match='active'):
        continuation.submit(config, True, '1')
    monkeypatch.setattr(continuation.subprocess, 'run', lambda argv, **kwargs: SimpleNamespace(
        stdout='', stderr='', returncode=0))
    ledger = root / 'fep/F64A/folded/w0_r0/attempts.json'
    ledger.parent.mkdir(parents=True)
    ledger.write_text(json.dumps([{}, {}]))
    with pytest.raises(ValueError, match='attempts'):
        continuation.submit(config, True, '1')


def test_failed_release_preserves_held_job_identity_and_prevents_duplicate(tmp_path, monkeypatch):
    config, root, _ = _staged_light(tmp_path, monkeypatch)

    def scheduler(argv, **kwargs):
        if argv[0] == 'qrls':
            raise subprocess.CalledProcessError(1, argv, stderr='synthetic release failure')
        return SimpleNamespace(stdout='789\n', stderr='', returncode=0)

    monkeypatch.setattr(continuation.subprocess, 'run', scheduler)
    with pytest.raises(subprocess.CalledProcessError):
        continuation.submit(config, True)
    assert json.loads((root.parent / 'first_light_submission.json').read_text())['job_id'] == '789'
    with pytest.raises(FileExistsError):
        continuation.submit(config, True)
