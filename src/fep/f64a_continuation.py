"""Isolated, resource-admitted 9→15 ns continuation of the completed F64A pilot."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from src.fep import f64a_extension as v1
from src.fep.pmx_engine import (
    GpuUnavailableError, ToolError, _GPU_BUSY_SIGNATURES, _run, dhdl_to_u_kn,
    gmx_command, mdrun_argv, verify_tools,
)
from src.prep.build import load_config

ROOT = v1.ROOT
DEFAULT_CONFIG = ROOT / "pilots/f64a_sampling_extension_v2/config.yaml"
BOUNDARY_POLICY = "freeze_v1_500_to_9500ps;carry_3500ps_exception;new_from_9501ps"


def repo_path(value: str) -> Path:
    """Resolve a configured path against the repository, never the caller's cwd."""
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def write_json(path: Path, value: dict) -> None:
    """Publish complete JSON atomically without overwriting an existing target."""
    serialized = json.dumps(value, indent=2, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(serialized)
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink()


def load_pilot(config: str | Path = DEFAULT_CONFIG) -> tuple[dict, dict]:
    """Validate the authorized arm and isolated roots; budgets live in admission."""
    path = Path(config).resolve()
    pilot = yaml.safe_load(path.read_text())
    expected = dict(pilot_id="f64a_sampling_extension_v2", status="authorized_2026-09-17",
                    variant="F64A", leg="folded", lambda_windows=20, replicates=3,
                    source_protocol="bf6841ccb3b9de79", baseline_ns=9., extension_ns=6.,
                    target_ns=15., baseline_total_time_ps=9500., energy_stride_ps=1.,
                    baseline_samples_per_window=9001, target_samples_per_window=15001,
                    first_light_extension_ps=10., trajectory_retention="estimates_only",
                    max_concurrent_tasks=8, gpu_type="L40S", walltime="12:00:00",
                    first_light_walltime="01:00:00", project="rise-batteries", cpu_slots=8,
                    prior_job_id=7589742)
    if any(pilot.get(k) != v for k, v in expected.items()):
        raise ValueError("v2 config differs from the authorized contract")
    cfg = load_config(repo_path(pilot["base_config"]))
    _, old = v1.load_pilot(repo_path(pilot["source_config"]))
    source, output = repo_path(pilot["source_production"]), repo_path(pilot["output_root"])
    if source != Path(old["_output_root"]) / "production":
        raise ValueError("source must be the configured completed v1 production root")
    if (output == source or output in source.parents or source in output.parents
            or output == ROOT / "results/fep" or ROOT / "results/fep" in output.parents
            or Path(old["_output_root"]) == output or Path(old["_output_root"]) in output.parents
            or ROOT / "results/pilots" not in output.parents):
        raise ValueError("v2 output must be an isolated pilot root")
    if (cfg["fep"]["lambda_windows"] != pilot["lambda_windows"]
            or cfg["fep"]["replicates"] != pilot["replicates"]
            or not cfg["fep"]["keep_disulfide_reduced"]
            or cfg["cluster"]["mdrun_ntomp"] != pilot["cpu_slots"]):
        raise ValueError("base protocol/state/resources differ from the selected arm")
    for key in ("runtime_headroom_multiplier", "storage_headroom_multiplier"):
        if not math.isfinite(pilot[key]) or pilot[key] < 1:
            raise ValueError("planning margins must be finite and at least one")
    pilot.update(_config_path=str(path), _output_root=str(output), _source_root=str(source),
                 _base_config_path=str(repo_path(pilot["base_config"])))
    return cfg, pilot


def code_identity(pilot: dict) -> dict:
    """Require committed simulation, analysis, environment and reference dependencies."""
    paths = [Path(__file__), ROOT / "src/analysis/f64a_continuation.py",
             ROOT / "src/analysis/f64a_blocks.py", ROOT / "src/analysis/f64a_sensitivity.py",
             ROOT / "src/fep/f64a_extension.py", ROOT / "src/fep/analyze.py",
             ROOT / "src/fep/pmx_engine.py", ROOT / "src/prep/build.py",
             ROOT / "scripts/scc_env.sh", ROOT / "workflow/Snakefile",
             Path(pilot["_config_path"]), Path(pilot["_base_config_path"]),
             *[repo_path(pilot[k]) for k in ("reference_primary", "reference_blocks",
                                           "reference_sensitivity", "source_config", "analysis_config")],
             *sorted((ROOT / "pilots/f64a_sampling_extension_v2").glob("*"))]
    if any(not p.is_file() for p in paths):
        raise FileNotFoundError("v2 dependency inventory is incomplete")
    names = sorted({str(p.resolve().relative_to(ROOT)) for p in paths})
    status = subprocess.run(["git", "status", "--porcelain", "--", *names], cwd=ROOT,
                            check=True, capture_output=True, text=True).stdout.strip()
    if status:
        raise RuntimeError(f"commit v2 dependencies before execution: {status}")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                            capture_output=True, text=True).stdout.strip()
    return {"git_commit": commit, "source_sha256": {p: v1.sha256_file(ROOT / p) for p in names}}


def source_files(folded: Path, window: int, rep: int) -> list[tuple[Path, Path]]:
    """Inventory the actual v1 working outputs, not absent original MDP/TPR files."""
    name = v1.task_name(window, rep)
    return [(folded / f"{name}.npz", Path(f"{name}.npz"))] + [
        (folded / name / p, Path(name) / p)
        for p in ("extension.tpr", "prod.cpt", "prod.log", "prod.edr", "dhdl.xvg",
                  "extension_protocol.sha")]


def scalar(text: str, name: str) -> float:
    """Read a unique scalar from GROMACS dump, refusing ambiguous output."""
    values = re.findall(r"^\s*" + re.escape(name) + r"\s*=\s*(\S+)", text, re.M)
    if len(values) != 1:
        raise ValueError(f"GROMACS dump needs one {name}, found {len(values)}")
    value = float(values[0])
    if not math.isfinite(value):
        raise ValueError(f"non-finite GROMACS scalar {name}")
    return value


def dump_tpr(cfg: dict, path: Path) -> str:
    """Capture stdout only so banners/command paths cannot contaminate model hashes."""
    result = subprocess.run([gmx_command(cfg), "dump", "-s", str(path)], cwd=path.parent,
                            check=True, capture_output=True, text=True)
    return result.stdout


def tpr_model_hash(dump: str) -> str:
    """Hash the complete dumped model/state with only duration (nsteps) excluded."""
    start = dump.find("inputrec:")
    if start < 0:
        raise ValueError("TPR dump has no inputrec section")
    lines = [line.rstrip() for line in dump[start:].splitlines()
             if not re.match(r"^\s*nsteps\s*=", line)]
    return hashlib.sha256('\n'.join(lines).encode()).hexdigest()


def check_tpr(cfg: dict, path: Path, window: int, end_ps: float,
              expected_model_hash: str | None = None) -> dict:
    """Verify the actual continuation TPR duration, lambda and output cadence."""
    dump = dump_tpr(cfg, path)
    values = {key: scalar(dump, key) for key in
              ("dt", "nsteps", "init-step", "tinit", "init-lambda-state", "nstdhdl",
               "nstxout", "nstvout", "nstfout", "nstxout-compressed")}
    end = values["tinit"] + (values["init-step"] + values["nsteps"]) * values["dt"]
    if (not math.isclose(end, end_ps, abs_tol=1e-6) or values["init-lambda-state"] != window
            or values["dt"] * values["nstdhdl"] != 1.
            or any(values[k] != 0 for k in ("nstxout", "nstvout", "nstfout", "nstxout-compressed"))):
        raise ValueError(f"{path}: incompatible TPR time/lambda/retention: {values}")
    model_hash = tpr_model_hash(dump)
    if expected_model_hash and model_hash != expected_model_hash:
        raise ValueError("TPR changed more than duration; model/state mismatch")
    values["model_sha256_excluding_nsteps"] = model_hash
    return values


def checkpoint_time(cfg: dict, path: Path, end_ps: float) -> None:
    """Verify checkpoint time independently of XVG's final row."""
    dump = _run([gmx_command(cfg), "dump", "-cp", str(path)], cwd=path.parent)
    if not math.isclose(scalar(dump, "t"), end_ps, abs_tol=1e-6):
        raise ValueError(f"{path}: checkpoint is not at {end_ps:g} ps")


def reconcile(xvg: Path, source_npz: Path, cfg: dict, pilot: dict,
              end_ps: float) -> tuple[np.ndarray, dict]:
    """Freeze the entire source matrix; permit only the two documented boundaries."""
    times = v1._xvg_times(xvg)
    if not np.array_equal(times, np.arange(0., end_ps + .5, pilot["energy_stride_ps"])):
        raise ValueError("XVG grid has gaps, repeats, wrong end time or cadence")
    parsed = dhdl_to_u_kn(xvg, cfg["fep"]["temperature_K"], pilot["lambda_windows"])
    if parsed.shape != (pilot["lambda_windows"], len(times)) or not np.isfinite(parsed).all():
        raise ValueError("invalid parsed XVG energies")
    with np.load(source_npz) as source:
        frozen = np.asarray(source["u_kn_window"])
        old_difference = float(source["append_boundary_max_abs_reduced_potential_difference"])
    historical = parsed[:, (times >= 500.) & (times <= pilot["baseline_total_time_ps"])]
    if historical.shape != frozen.shape:
        raise ValueError("source history shape mismatch")
    old_column = int(3500. - 500.)
    measured_old = float(np.max(abs(historical[:, old_column] - frozen[:, old_column])))
    if measured_old != old_difference:
        raise ValueError("3500 ps boundary does not match recorded v1 exception")
    exact = np.ones(frozen.shape[1], dtype=bool)
    exact[old_column] = False
    continued = end_ps > pilot["baseline_total_time_ps"]
    if continued:
        exact[-1] = False
    if not np.array_equal(historical[:, exact], frozen[:, exact]):
        raise ValueError("historical energies changed outside documented boundaries")
    new_difference = float(np.max(abs(historical[:, -1] - frozen[:, -1]))) if continued else 0.
    new = parsed[:, times > pilot["baseline_total_time_ps"]]
    if new.shape[1] != int(end_ps - pilot["baseline_total_time_ps"]):
        raise ValueError("wrong extension sample count")
    return new, {"3500": measured_old, "9500": new_difference}


def validate_source(config: str | Path = DEFAULT_CONFIG) -> dict:
    """Read-only validation and checksumming of all source files and ancestry."""
    cfg, pilot = load_pilot(config)
    root = Path(pilot["_source_root"])
    folded = root / "fep/F64A/folded"
    _, old = v1.load_pilot(repo_path(pilot["source_config"]))
    old_stage = v1._load_stage_manifest(root, old, require_current_code=False)
    sensitivity = json.loads(repo_path(pilot["reference_sensitivity"]).read_text())
    old_manifest_hash, = [h for p, h in sensitivity["input_sha256"].items()
                          if Path(p).name == "stage_manifest.json"]
    if (v1.sha256_file(root / "stage_manifest.json") != old_manifest_hash
            or old_stage["code_identity"] != sensitivity["simulation_code_identity"]):
        raise ValueError("v1 stage manifest differs from frozen ancestry")
    expected_hashes = {Path(p).name: h for p, h in sensitivity["input_sha256"].items()
                       if p.endswith(".npz")}
    entries, tpr_records = [], []
    for window in range(pilot["lambda_windows"]):
        for rep in range(pilot["replicates"]):
            name = v1.task_name(window, rep)
            for source, relative in source_files(folded, window, rep):
                if not source.is_file() or source.stat().st_size == 0:
                    raise FileNotFoundError(f"v2 source missing/empty: {source}")
                entries.append({"source": str(source), "relative": str(relative),
                                "sha256": v1.sha256_file(source), "size_bytes": source.stat().st_size})
            npz = folded / f"{name}.npz"
            if v1.sha256_file(npz) != expected_hashes.get(npz.name):
                raise ValueError(f"{name}: NPZ differs from frozen sensitivity input")
            v1._validate_output_npz(npz, old, old["extension_ns"] * 1000., old_stage)
            if (folded / name / "extension_protocol.sha").read_text().strip() != pilot["source_protocol"]:
                raise ValueError("v1 continuation stamp mismatch")
            reconcile(folded / name / "dhdl.xvg", npz, cfg, pilot, pilot["baseline_total_time_ps"])
            original = root / "source/F64A/folded" / name / "prod.tpr"
            original_entry, = [e for e in old_stage["files"] if e["path"] == str(original.relative_to(root))]
            if v1.sha256_file(original) != original_entry["sha256"]:
                raise ValueError("original v1 TPR changed since v1 staging")
            model_hash = tpr_model_hash(dump_tpr(cfg, original))
            tpr_records.append({"task": name, "values": check_tpr(
                cfg, folded / name / "extension.tpr", window, pilot["baseline_total_time_ps"], model_hash)})
            checkpoint_time(cfg, folded / name / "prod.cpt", pilot["baseline_total_time_ps"])
    for key, filename in (("reference_primary", "f64a_extension.json"),
                          ("reference_blocks", "f64a_blocks_v1.json"),
                          ("reference_sensitivity", "f64a_selection_sensitivity_v1.json")):
        cluster = root / "analysis" / filename
        if json.loads(cluster.read_text()) != json.loads(repo_path(pilot[key]).read_text()):
            raise ValueError(f"cluster report differs from frozen {filename}")
        entries.append({"source": str(cluster), "relative": "ancestry/" + filename,
                        "sha256": v1.sha256_file(cluster), "size_bytes": cluster.stat().st_size})
    manifest = root / "stage_manifest.json"
    entries.append({"source": str(manifest), "relative": "ancestry/v1_stage_manifest.json",
                    "sha256": v1.sha256_file(manifest), "size_bytes": manifest.stat().st_size})
    if any(v1.sha256_file(e["source"]) != e["sha256"] for e in entries):
        raise ValueError("v1 source changed during validation")
    return {"pilot_id": pilot["pilot_id"], "files": entries, "tasks": 60,
            "bytes": sum(e["size_bytes"] for e in entries), "tpr_checks": tpr_records,
            "source_code_identity": old_stage["code_identity"]}


def parse_accounting(text: str, job_id: int) -> dict:
    """Require all 60 successful array tasks, not a qstat elapsed-time estimate."""
    records = []
    for block in re.split(r"^=+\s*$", text, flags=re.M):
        fields = dict(re.findall(r"^(\S+)\s+(.+?)\s*$", block, re.M))
        if "taskid" not in fields:
            continue
        if int(fields["jobnumber"]) != job_id:
            raise ValueError("accounting belongs to another job")
        if int(fields["failed"]) != 0 or int(fields["exit_status"]) != 0:
            raise ValueError("prior array has failed accounting attempts; review before admission")
        task, wall = int(fields["taskid"]), float(fields["ru_wallclock"])
        if not math.isfinite(wall) or wall <= 0:
            raise ValueError("invalid accounted wallclock")
        cpu = float(fields["cpu"])
        if not math.isfinite(cpu) or cpu < 0:
            raise ValueError("invalid accounted CPU usage")
        records.append((task, wall, cpu))
    if len(records) != 60 or sorted(r[0] for r in records) != list(range(1, 61)):
        raise ValueError("accounting must contain exactly 60 successful tasks")
    return {"gpu_hours": sum(r[1] for r in records) / 3600.,
            "maximum_task_seconds": max(r[1] for r in records),
            "cpu_hours": sum(r[2] for r in records) / 3600.}


def admit(config: str | Path, qacct_file: Path, quota_file: Path, balance_file: Path,
          available_gib: float, su_balance: float, su_rate: float, gpu_hour_cap: float,
          su_cap: float, max_attempts: int) -> Path:
    """Require explicit resource caps and measured evidence before any staging/GPU use."""
    _, pilot = load_pilot(config)
    identity = code_identity(pilot)
    values = (available_gib, su_balance, su_rate, gpu_hour_cap, su_cap)
    if any(not math.isfinite(v) or v <= 0 for v in values) or type(max_attempts) is not int or max_attempts < 1:
        raise ValueError("resource values must be positive and finite; attempts must be an integer")
    evidence = {}
    for label, path in (("qacct", qacct_file), ("quota", quota_file), ("balance", balance_file)):
        if not path.is_file() or not path.read_text().strip():
            raise ValueError(f"missing resource evidence: {label}")
        evidence[label] = {"path": str(path.resolve()), "sha256": v1.sha256_file(path),
                           "text": path.read_text()}
    accounting = parse_accounting(evidence["qacct"]["text"], pilot["prior_job_id"])
    source = validate_source(config)
    # Full frozen staging + private growing append sets + output matrices, plus a
    # conservative margin; first-light storage is included by doubling this total.
    ratio = pilot["target_ns"] / pilot["baseline_ns"]
    storage = math.ceil(pilot["storage_headroom_multiplier"] * (
        source["bytes"] * (1 + ratio) + 60 * pilot["lambda_windows"] * pilot["target_samples_per_window"] * 8))
    first_seconds = 3600
    task_seconds = min(12 * 3600, math.floor((gpu_hour_cap * 3600 - max_attempts * first_seconds)
                                           / (60 * max_attempts)))
    if task_seconds < accounting["maximum_task_seconds"] * pilot["runtime_headroom_multiplier"]:
        raise ValueError("GPU-hour cap cannot accommodate measured task runtime plus declared headroom")
    reserved_hours = max_attempts * (first_seconds + 60 * task_seconds) / 3600.
    reserved_su = reserved_hours * pilot["cpu_slots"] * su_rate
    if storage > available_gib * 1024**3 or reserved_su > min(su_cap, su_balance):
        raise ValueError("insufficient quota headroom or approved/project SU budget")
    output = Path(pilot["_output_root"]) / "admission.json"
    write_json(output, {"schema_version": 1, "pilot_id": pilot["pilot_id"],
                       "created_utc": datetime.now(timezone.utc).isoformat(),
                       "code_identity": identity, "pilot_config_sha256": v1.sha256_file(config),
                       "source_inventory": source, "evidence": evidence, "accounting": accounting,
                       "available_gib": available_gib, "su_balance": su_balance, "su_rate": su_rate,
                       "gpu_hour_cap": gpu_hour_cap, "su_cap": su_cap, "max_attempts": max_attempts,
                       "first_light_seconds": first_seconds, "task_seconds": task_seconds,
                       "maximum_reserved_gpu_hours": reserved_hours, "maximum_reserved_su": reserved_su,
                       "required_storage_bytes": storage})
    print(f"Admitted resources: {output}; GPU cap {reserved_hours:.2f} h, SU cap {reserved_su:.2f}")
    return output


def admission(pilot: dict, current: bool = True) -> dict:
    """Load the immutable resource approval and enforce its code/config identity."""
    value = json.loads((Path(pilot["_output_root"]) / "admission.json").read_text())
    if value["pilot_id"] != pilot["pilot_id"] or value["pilot_config_sha256"] != v1.sha256_file(pilot["_config_path"]):
        raise ValueError("admission pilot/config mismatch")
    if current and value["code_identity"] != code_identity(pilot):
        raise ValueError("code changed after resource admission")
    return value


def stage(config: str | Path = DEFAULT_CONFIG, first_light: bool = False) -> Path:
    """Stage checksum-frozen full append files in a new namespace; no overwrites."""
    _, pilot = load_pilot(config)
    approved = admission(pilot)
    mode = "first_light" if first_light else "production"
    if not first_light:
        check_first_light(config)
    root = Path(pilot["_output_root"]) / mode
    if root.exists():
        load_stage(pilot, first_light)
        return root
    files = approved["source_inventory"]["files"]
    if first_light:
        files = [e for e in files if e["relative"] == "w0_r0.npz" or e["relative"].startswith(("w0_r0/", "ancestry/"))]
    if any(v1.sha256_file(e["source"]) != e["sha256"] for e in files):
        raise ValueError("source changed since resource admission")
    entries = []
    for entry in files:
        destination = root / "source/F64A/folded" / entry["relative"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        v1._copy_atomic(Path(entry["source"]), destination)
        if v1.sha256_file(destination) != entry["sha256"]:
            raise ValueError("copied source checksum mismatch")
        entries.append({"path": str(destination.relative_to(root)), **entry})
    write_json(root / "stage_manifest.json", {"pilot_id": pilot["pilot_id"], "mode": mode,
               "code_identity": approved["code_identity"], "admission_sha256":
               v1.sha256_file(Path(pilot["_output_root"]) / "admission.json"), "files": entries,
               "source_code_identity": approved["source_inventory"]["source_code_identity"]})
    print(f"Staged {mode}: {root}")
    return root


def load_stage(pilot: dict, first_light: bool, current: bool = True,
               task: str | None = None) -> tuple[Path, dict]:
    """Check admission, stage identity and hashes, optionally for one task only."""
    mode = "first_light" if first_light else "production"
    root = Path(pilot["_output_root"]) / mode
    manifest = json.loads((root / "stage_manifest.json").read_text())
    approved = admission(pilot, current)
    if (manifest["pilot_id"] != pilot["pilot_id"] or manifest["mode"] != mode
            or manifest["code_identity"] != approved["code_identity"]
            or manifest["admission_sha256"] != v1.sha256_file(Path(pilot["_output_root"]) / "admission.json")):
        raise ValueError("stage identity/admission mismatch")
    for entry in manifest["files"]:
        if task and entry["relative"] != task + ".npz" and not entry["relative"].startswith(task + "/"):
            continue
        if v1.sha256_file(root / entry["path"]) != entry["sha256"]:
            raise ValueError("frozen staged source changed")
    return root, manifest


def validate_output(path: Path, pilot: dict, manifest: dict, extension_ps: float,
                    source_name: str | None = None) -> None:
    """Validate complete metadata, exact ancestry columns and recorded boundaries."""
    logical = Path(source_name or path.name)
    staged = path.parents[3] / "source/F64A/folded" / logical.name
    with np.load(staged) as source, np.load(path) as data:
        expected_shape = (pilot["lambda_windows"], pilot["baseline_samples_per_window"] + int(extension_ps))
        u = data["u_kn_window"]
        if u.shape != expected_shape or not np.isfinite(u).all():
            raise ValueError("invalid output energy shape/values")
        if not np.array_equal(u[:, :pilot["baseline_samples_per_window"]], source["u_kn_window"]):
            raise ValueError("output changed frozen 9 ns source columns")
        expected = {"pilot_id": pilot["pilot_id"], "provenance": "gromacs_pmx", "n_states": 20,
                    "lambda_index": int(logical.stem.split('_')[0][1:]),
                    "replicate": int(logical.stem.split('_r')[1]),
                    "protocol": v1._extension_protocol(pilot, extension_ps),
                    "source_protocol": pilot["source_protocol"], "source_npz_sha256": v1.sha256_file(staged),
                    "code_identity_json": json.dumps(manifest["code_identity"], sort_keys=True),
                    "source_code_identity_json": json.dumps(manifest["source_code_identity"], sort_keys=True),
                    "admission_sha256": manifest["admission_sha256"], "append_boundary_policy": BOUNDARY_POLICY}
        for key, value in expected.items():
            if str(data[key]) != str(value):
                raise ValueError(f"output metadata mismatch: {key}")
        deltas = json.loads(str(data["append_boundary_differences_json"]))
        if (deltas["3500"] != float(source["append_boundary_max_abs_reduced_potential_difference"])
                or not math.isfinite(deltas["9500"]) or deltas["9500"] < 0):
            raise ValueError("invalid boundary difference records")


def run_task(window: int, rep: int, config: str | Path = DEFAULT_CONFIG,
             first_light: bool = False) -> Path:
    """Continue one private append set, with a finite admission-backed attempt budget."""
    cfg, pilot = load_pilot(config)
    v1._validate_task(pilot, window, rep)
    if first_light and (window, rep) != (0, 0):
        raise ValueError("first light is w0/r0 only")
    if not first_light:
        check_first_light(config)
    name = v1.task_name(window, rep)
    root, manifest = load_stage(pilot, first_light, task=name)
    approved = admission(pilot)
    extension = pilot["first_light_extension_ps"] if first_light else pilot["extension_ns"] * 1000.
    source = root / "source/F64A/folded"
    out = root / "fep/F64A/folded" / (name + ".npz")
    if out.exists():
        validate_output(out, pilot, manifest, extension)
        return out
    if not os.environ.get("JOB_ID"):
        raise RuntimeError("GPU continuation requires an SGE job; use the guarded submit command")
    mode = "first_light" if first_light else "production"
    submissions = [json.loads(p.read_text()) for p in Path(pilot["_output_root"]).glob(mode + "_submission*.json")]
    if not any(s["job_id"].split(".")[0] == os.environ["JOB_ID"] for s in submissions):
        raise RuntimeError("job is not the recorded guarded submission; refusing unbudgeted duplicate")
    run = out.parent / name
    run.mkdir(parents=True, exist_ok=True)
    with (run / ".task.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        attempts_path = run / "attempts.json"
        attempts = json.loads(attempts_path.read_text()) if attempts_path.exists() else []
        if len(attempts) >= approved["max_attempts"]:
            raise ToolError("approved task attempts exhausted; no automatic rescheduling")
        attempts.append({"job_id": os.environ["JOB_ID"], "utc": datetime.now(timezone.utc).isoformat()})
        temporary = attempts_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(attempts))
        temporary.replace(attempts_path)
        protocol = v1._extension_protocol(pilot, extension)
        stamp = run / "extension_protocol.sha"
        if stamp.exists() and stamp.read_text().strip() != protocol:
            raise ValueError("run protocol changed; refusing resume")
        if not stamp.exists() and (run / "prod.cpt").exists():
            raise ValueError("unstamped checkpoint; refusing resume")
        stamp.write_text(protocol + "\n")
        for filename in ("prod.cpt", "prod.log", "prod.edr", "dhdl.xvg"):
            if not (run / filename).exists():
                v1._copy_atomic(source / name / filename, run / filename)
        tpr = run / "extension.tpr"
        end = pilot["baseline_total_time_ps"] + extension
        expected_model, = [entry["values"]["model_sha256_excluding_nsteps"]
                           for entry in approved["source_inventory"]["tpr_checks"] if entry["task"] == name]
        if not tpr.exists():
            temporary_tpr = run / f"extension.building.{os.getpid()}.tpr"
            _run([gmx_command(cfg), "convert-tpr", "-s", str(source / name / "extension.tpr"),
                  "-extend", str(extension), "-o", str(temporary_tpr)], cwd=run)
            check_tpr(cfg, temporary_tpr, window, end, expected_model)
            temporary_tpr.replace(tpr)
        else:
            check_tpr(cfg, tpr, window, end, expected_model)
        seconds = approved["first_light_seconds"] if first_light else approved["task_seconds"]
        argv = [*mdrun_argv(cfg), "-s", "extension.tpr", "-deffnm", "prod", "-cpi", "prod.cpt",
                "-append", "-cpo", "prod.cpt", "-cpt", str(cfg["fep"]["checkpoint_interval_min"]),
                "-dhdl", "dhdl.xvg", "-maxh", str(seconds / 3600.)]
        try:
            result = subprocess.run(argv, cwd=run, capture_output=True, text=True, timeout=seconds,
                                    env={**os.environ, "OMP_NUM_THREADS": str(pilot["cpu_slots"])})
        except subprocess.TimeoutExpired as exc:
            raise ToolError("admitted walltime exhausted; resume only within remaining attempt budget") from exc
        output = (result.stdout or "") + (result.stderr or "")
        print(output[-8000:], flush=True)
        if result.returncode:
            if any(s in output for s in _GPU_BUSY_SIGNATURES) and len(attempts) < approved["max_attempts"]:
                raise GpuUnavailableError("scheduler assigned an unusable GPU; finite-budget reschedule")
            raise ToolError(f"mdrun failed ({result.returncode}): {output[-4000:]}")
        checkpoint_time(cfg, run / "prod.cpt", end)
        new, deltas = reconcile(run / "dhdl.xvg", source / (name + ".npz"), cfg, pilot, end)
        with np.load(source / (name + ".npz")) as baseline:
            combined = np.concatenate([baseline["u_kn_window"], new], axis=1)
        tmp = out.with_suffix(".tmp.npz")
        np.savez(tmp, u_kn_window=combined, lambda_index=window, replicate=rep, n_states=20,
                 provenance="gromacs_pmx", pilot_id=pilot["pilot_id"], protocol=protocol,
                 source_protocol=pilot["source_protocol"], source_npz_sha256=v1.sha256_file(source / (name + ".npz")),
                 code_identity_json=json.dumps(manifest["code_identity"], sort_keys=True),
                 source_code_identity_json=json.dumps(manifest["source_code_identity"], sort_keys=True),
                 admission_sha256=manifest["admission_sha256"], append_boundary_policy=BOUNDARY_POLICY,
                 append_boundary_differences_json=json.dumps(deltas, sort_keys=True))
        validate_output(tmp, pilot, manifest, extension, source_name=out.name)
        tmp.replace(out)
    print(f"Wrote {out}: {combined.shape}")
    return out


def check_first_light(config: str | Path = DEFAULT_CONFIG) -> Path:
    """Require a complete, finite, metadata-verified real first-light result."""
    _, pilot = load_pilot(config)
    root, manifest = load_stage(pilot, True, task="w0_r0")
    out = root / "fep/F64A/folded/w0_r0.npz"
    validate_output(out, pilot, manifest, pilot["first_light_extension_ps"])
    print(f"First light passed: {out} (20, 9011)")
    return out


def validate_outputs(config: str | Path = DEFAULT_CONFIG) -> tuple[Path, dict]:
    """Validate all 60 immutable outputs without requiring old simulation Git HEAD."""
    _, pilot = load_pilot(config)
    root, manifest = load_stage(pilot, False, current=False)
    fep = root / "fep"
    for window in range(pilot["lambda_windows"]):
        for rep in range(pilot["replicates"]):
            validate_output(fep / "F64A/folded" / (v1.task_name(window, rep) + ".npz"),
                            pilot, manifest, pilot["extension_ns"] * 1000.)
    print("Validated 60 v2 outputs: (20, 15001), exact frozen 9 ns prefixes")
    return fep, manifest


def submit(config: str | Path, first_light: bool, retry_tasks: str | None = None) -> str:
    """Submit only staged, resource-admitted jobs; record IDs and reject duplicate launch."""
    _, pilot = load_pilot(config)
    approved = admission(pilot)
    mode = "first_light" if first_light else "production"
    load_stage(pilot, first_light)
    if not first_light:
        check_first_light(config)
    marker = Path(pilot["_output_root"]) / (mode + "_submission.json")
    seconds = approved["first_light_seconds"] if first_light else approved["task_seconds"]
    walltime = f"{seconds // 3600:02d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"
    script = ROOT / "pilots/f64a_sampling_extension_v2" / ("submit_first_light.sh" if first_light else "submit_array.sh")
    (ROOT / "logs/f64a_continuation").mkdir(parents=True, exist_ok=True)
    output_root = Path(pilot["_output_root"])
    with (output_root / ".submission.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        extra, tasks = [], [1] if first_light else list(range(1, 61))
        previous = list(output_root.glob(mode + "_submission*.json"))
        if not retry_tasks and previous:
            raise FileExistsError("already submitted; inspect recorded job, do not duplicate it")
        if retry_tasks:
            if not previous:
                raise ValueError("retry requires a recorded prior submission")
            tasks = [int(value) for value in retry_tasks.split(',')]
            if (len(set(tasks)) != len(tasks) or not tasks or any(t < 1 or t > 60 for t in tasks)
                    or first_light and tasks != [1]):
                raise ValueError("invalid retry task IDs")
            if sorted(tasks) != list(range(min(tasks), max(tasks) + 1)):
                raise ValueError("SCC retry task IDs must be contiguous; submit separate retries after each finishes")
            jobs = subprocess.run(["qstat", "-u", os.environ["USER"]], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout
            active = set(re.findall(r"^\s*(\d+)\s", jobs, re.M))
            if any(json.loads(p.read_text())["job_id"].split('.')[0] in active for p in previous):
                raise RuntimeError("prior submission still active; do not duplicate tasks")
            root = output_root / mode / "fep/F64A/folded"
            for task in tasks:
                name = v1.task_name((task - 1) // pilot['replicates'], (task - 1) % pilot['replicates'])
                if (root / (name + '.npz')).exists():
                    raise ValueError("retry includes an existing result; diagnose rather than overwrite")
                ledger = root / name / 'attempts.json'
                if ledger.exists() and len(json.loads(ledger.read_text())) >= approved['max_attempts']:
                    raise ValueError("retry exceeds admitted attempts")
            marker = output_root / f"{mode}_submission_retry_{len(previous)}.json"
            if not first_light:
                extra = ['-t', f'{min(tasks)}-{max(tasks)}']
        # Hold until provenance exists: a promptly dispatched task must never race
        # the submission-record write and fail its guarded JOB_ID check.
        result = subprocess.run(["qsub", "-h", "-terse", "-l", "h_rt=" + walltime, *extra, str(script)],
                                cwd=ROOT, check=True, capture_output=True, text=True)
        job = result.stdout.strip()
        if not re.fullmatch(r"\d+(?:\.[0-9,:-]+)?", job):
            raise RuntimeError(f"qsub returned unexpected job identity; inspect before retry: {job}")
        print(f"qsub accepted {job}; recording submission", flush=True)
        write_json(marker, {"job_id": job, "walltime": walltime, "mode": mode, "tasks": tasks,
                            "admission_sha256": v1.sha256_file(output_root / "admission.json")})
        subprocess.run(["qrls", job.split('.')[0]], cwd=ROOT, check=True,
                       capture_output=True, text=True)
    print(f"Submitted {mode}: {job}; h_rt={walltime}")
    return job


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate-source", "verify-tools", "check-first-light", "validate-outputs"):
        sub.add_parser(name)
    for name in ("stage", "submit"):
        command = sub.add_parser(name)
        command.add_argument("--first-light", action="store_true")
        if name == "submit":
            command.add_argument("--retry-tasks", help="explicit failed task IDs after prior job finishes")
    run = sub.add_parser("run")
    run.add_argument("--window", type=int, required=True)
    run.add_argument("--rep", type=int, required=True)
    run.add_argument("--first-light", action="store_true")
    approve = sub.add_parser("admit")
    for name in ("qacct-file", "quota-file", "balance-file"):
        approve.add_argument("--" + name, type=Path, required=True)
    for name in ("available-gib", "su-balance", "su-rate", "gpu-hour-cap", "su-cap"):
        approve.add_argument("--" + name, type=float, required=True)
    approve.add_argument("--max-attempts", type=int, required=True)
    args = parser.parse_args()
    try:
        if args.command == "validate-source":
            value = validate_source(args.config)
            print(f"Validated {value['tasks']} tasks / {len(value['files'])} files / {value['bytes']} bytes")
        elif args.command == "verify-tools":
            cfg, _ = load_pilot(args.config)
            print(json.dumps(verify_tools(cfg), indent=2))
        elif args.command == "admit":
            admit(args.config, args.qacct_file, args.quota_file, args.balance_file, args.available_gib,
                  args.su_balance, args.su_rate, args.gpu_hour_cap, args.su_cap, args.max_attempts)
        elif args.command == "stage":
            stage(args.config, args.first_light)
        elif args.command == "submit":
            submit(args.config, args.first_light, args.retry_tasks)
        elif args.command == "run":
            run_task(args.window, args.rep, args.config, args.first_light)
        elif args.command == "check-first-light":
            check_first_light(args.config)
        else:
            validate_outputs(args.config)
    except GpuUnavailableError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(99)


if __name__ == "__main__":
    main()
