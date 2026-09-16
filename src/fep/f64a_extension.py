"""Isolated F64A continuation pilot: extend the retained 3 ns paths to 9 ns.

This module deliberately does not use :func:`src.fep.window.run_window`.  A normal
window starts from the minimised system and writes below ``results/fep``; this pilot
continues the archived checkpoint and writes below a separate, preregistered root.
The historical archive is read only during staging and is never used as an mdrun
working directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from src.fep.analyze import _KB_KCAL, _logmeanexp, decorrelate_window, solve_leg_mbar
from src.fep.pmx_engine import (
    GpuUnavailableError,
    ToolError,
    _GPU_BUSY_SIGNATURES,
    _run,
    dhdl_to_u_kn,
    gmx_command,
    mdrun_argv,
    verify_tools,
)
from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PILOT_CONFIG = ROOT / "pilots" / "f64a_sampling_extension_v1" / "config.yaml"
PILOT_GIT_PATHS = (
    "src/fep/f64a_extension.py",
    "pilots/f64a_sampling_extension_v1",
    "workflow/Snakefile",
    "config/pipeline.yaml",
)


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _copy_atomic(source: Path, destination: Path) -> None:
    """Copy through a sibling temporary file so a killed job leaves no partial target."""
    temporary = destination.with_name(destination.name + ".copying")
    shutil.copy2(source, temporary)
    temporary.replace(destination)


def code_identity() -> dict:
    """Return the repository commit and exact pilot-module hash."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--", *PILOT_GIT_PATHS], cwd=ROOT,
            check=True, capture_output=True, text=True,
        ).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError("pilot execution requires an identifiable Git checkout") from exc
    return {
        "git_commit": commit,
        "pilot_code_sha256": sha256_file(Path(__file__).resolve()),
        "pilot_paths_clean": not status,
        "pilot_git_status": status,
    }


def _require_committed_code() -> dict:
    identity = code_identity()
    if not identity["pilot_paths_clean"]:
        raise RuntimeError(
            "pilot code/config is modified or untracked; commit it before staging: "
            f"{identity['pilot_git_status']}"
        )
    return identity


def adjacent_discrepancies(windows: list[np.ndarray], kt: float) -> np.ndarray:
    """Signed forward-minus-reverse adjacent discrepancies in kcal/mol."""
    return np.asarray([
        (-_logmeanexp(-(windows[k][k + 1] - windows[k][k]))
         - _logmeanexp(windows[k + 1][k + 1] - windows[k + 1][k])) * kt
        for k in range(len(windows) - 1)
    ])


def load_pilot(path: str | Path = DEFAULT_PILOT_CONFIG) -> tuple[dict, dict]:
    """Load the base pipeline and the fixed pilot settings, then validate the arm."""
    path = Path(path).resolve()
    pilot = yaml.safe_load(path.read_text())
    base = Path(pilot["base_config"])
    if not base.is_absolute():
        base = ROOT / base
    cfg = load_config(base)

    required = {
        "pilot_id": "f64a_sampling_extension_v1",
        "status": "authorized_2026-09-15",
        "variant": "F64A",
        "leg": "folded",
        "lambda_windows": 20,
        "replicates": 3,
        "baseline_ns": 3.0,
        "extension_ns": 6.0,
        "target_ns": 9.0,
        "baseline_samples_per_window": 3001,
        "target_samples_per_window": 9001,
        "source_protocol": "822108e9db71124d",
        "baseline_total_time_ps": 3500.0,
        "energy_stride_ps": 1.0,
        "first_light_extension_ps": 10.0,
        "trajectory_retention": "estimates_only",
        "max_concurrent_tasks": 8,
        "gpu_type": "L40S",
        "walltime": "12:00:00",
        "project": "rise-batteries",
    }
    wrong = {key: (pilot.get(key), want) for key, want in required.items()
             if pilot.get(key) != want}
    if wrong:
        raise ValueError(f"pilot config differs from its registered contract: {wrong}")
    if int(cfg["fep"]["lambda_windows"]) != pilot["lambda_windows"]:
        raise ValueError("base config lambda count no longer matches the registered pilot")
    if int(cfg["fep"]["replicates"]) != pilot["replicates"]:
        raise ValueError("base config replicate count no longer matches the registered pilot")

    output_root = (ROOT / pilot["output_root"]).resolve()
    historical = (ROOT / "results" / "fep").resolve()
    if output_root == historical or historical in output_root.parents:
        raise ValueError(f"pilot output {output_root} overlaps historical results {historical}")
    pilot["_config_path"] = str(path)
    pilot["_base_config_path"] = str(base.resolve())
    pilot["_output_root"] = str(output_root)
    return cfg, pilot


def task_name(window: int, rep: int) -> str:
    return f"w{window}_r{rep}"


def _validate_task(pilot: dict, window: int, rep: int) -> None:
    if not 0 <= window < int(pilot["lambda_windows"]):
        raise ValueError(f"window {window} is outside 0..{pilot['lambda_windows'] - 1}")
    if not 0 <= rep < int(pilot["replicates"]):
        raise ValueError(f"replicate {rep} is outside 0..{pilot['replicates'] - 1}")


def _source_files(source_folded: Path, window: int, rep: int) -> list[tuple[Path, Path]]:
    name = task_name(window, rep)
    run = source_folded / name
    return [
        (source_folded / f"{name}.npz", Path(f"{name}.npz")),
        (run / "prod.tpr", Path(name) / "prod.tpr"),
        (run / "prod.cpt", Path(name) / "prod.cpt"),
        (run / "prod.mdp", Path(name) / "prod.mdp"),
        # GROMACS only appends from a checkpoint when the prior appendable files match
        # the checksums stored in it. Copy them into the isolated working directory;
        # otherwise -append correctly refuses the continuation.
        (run / "prod.log", Path(name) / "prod.log"),
        (run / "prod.edr", Path(name) / "prod.edr"),
        (run / "dhdl.xvg", Path(name) / "dhdl.xvg"),
        (run / "protocol.sha", Path(name) / "protocol.sha"),
    ]


def _validate_baseline_npz(path: Path, pilot: dict, window: int) -> None:
    with np.load(path) as data:
        u = np.asarray(data["u_kn_window"])
        if u.shape != (pilot["lambda_windows"], pilot["baseline_samples_per_window"]):
            raise ValueError(f"{path}: baseline shape {u.shape} is not (20, 3001)")
        if int(data["lambda_index"]) != window or int(data["n_states"]) != 20:
            raise ValueError(f"{path}: lambda metadata does not identify window {window}/20")
        if str(data["provenance"]) != "gromacs_pmx":
            raise ValueError(f"{path}: baseline provenance is not gromacs_pmx")
        if str(data["protocol"]) != pilot["source_protocol"]:
            raise ValueError(f"{path}: baseline protocol is not {pilot['source_protocol']}")
        if not np.isfinite(u).all():
            raise ValueError(f"{path}: baseline contains non-finite energies")


def _mdp_values(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        body = line.split(";", 1)[0].strip()
        if "=" in body:
            key, value = body.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def _validate_baseline_mdp(path: Path, window: int) -> None:
    values = _mdp_values(path)
    expected = {
        "dt": "0.002",
        "nsteps": "1750000",
        "nstdhdl": "500",
        "nstxout-compressed": "0",
        "init-lambda-state": str(window),
    }
    wrong = {key: (values.get(key), want) for key, want in expected.items()
             if values.get(key) != want}
    if wrong:
        raise ValueError(f"{path}: baseline MDP does not match the continuation contract: {wrong}")


def _xvg_times(path: Path) -> np.ndarray:
    values = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        values.append(float(stripped.split()[0]))
    return np.asarray(values)


def reconcile_energy_history(xvg: Path, baseline_npz: Path, cfg: dict, pilot: dict,
                             expected_end_ps: float) -> np.ndarray:
    """Validate the archived 0.5–3.5 ns energies and return only new columns.

    The exact integer-picosecond time grid detects gaps, duplicate restart records, and
    truncated appends. Re-parsing the archived 500–3500 ps records and comparing them
    exactly with the frozen NPZ binds the XVG/checkpoint append set to that NPZ.
    """
    stride = float(pilot["energy_stride_ps"])
    expected_times = np.arange(0.0, expected_end_ps + 0.5 * stride, stride)
    times = _xvg_times(xvg)
    if not np.array_equal(times, expected_times):
        mismatch = next(
            (i for i, (got, want) in enumerate(zip(times, expected_times)) if got != want),
            min(len(times), len(expected_times)),
        )
        got = times[mismatch] if mismatch < len(times) else "<missing>"
        want = expected_times[mismatch] if mismatch < len(expected_times) else "<none>"
        raise ValueError(
            f"{xvg}: time grid is not exact 0..{expected_end_ps:g} ps by {stride:g} ps; "
            f"first mismatch at row {mismatch}: got {got}, expected {want}"
        )

    parsed = dhdl_to_u_kn(xvg, float(cfg["fep"]["temperature_K"]), 20)
    if parsed.shape[1] != len(times):
        raise ValueError(f"{xvg}: {len(times)} times but {parsed.shape[1]} energy rows")
    baseline_start_ps = (
        float(pilot["baseline_total_time_ps"]) - float(pilot["baseline_ns"]) * 1000.0
    )
    baseline_mask = ((times >= baseline_start_ps) &
                     (times <= float(pilot["baseline_total_time_ps"])))
    with np.load(baseline_npz) as baseline:
        frozen = np.asarray(baseline["u_kn_window"])
    reparsed = parsed[:, baseline_mask]
    if reparsed.shape != frozen.shape or not np.array_equal(reparsed, frozen):
        detail = f"shape {reparsed.shape} vs {frozen.shape}"
        if reparsed.shape == frozen.shape:
            detail = f"maximum absolute difference {np.max(np.abs(reparsed - frozen))}"
        raise ValueError(
            f"{xvg}: archived 500–3500 ps energies do not exactly match {baseline_npz} "
            f"({detail})"
        )

    new_mask = times > float(pilot["baseline_total_time_ps"])
    new = parsed[:, new_mask]
    expected_new = int(round(
        (expected_end_ps - float(pilot["baseline_total_time_ps"])) / stride
    ))
    if new.shape != (20, expected_new):
        raise ValueError(f"{xvg}: new energy shape {new.shape}, expected (20, {expected_new})")
    return new


def validate_source(source_folded: str | Path,
                    config: str | Path = DEFAULT_PILOT_CONFIG,
                    first_light: bool = False) -> dict:
    """Read-only validation of the archived inputs required by this arm."""
    cfg, pilot = load_pilot(config)
    source_folded = Path(source_folded).resolve()
    tasks = [(0, 0)] if first_light else [
        (w, r) for w in range(pilot["lambda_windows"])
        for r in range(pilot["replicates"])
    ]
    frozen_manifest = json.loads((ROOT / "data/forensic_archive_manifest.json").read_text())
    known_npz = {
        row["path"]: row["sha256"] for row in frozen_manifest["files"]
        if row["root"] == "main_archive" and row["path"].startswith("fep/F64A/folded/")
    }
    total_bytes = 0
    for window, rep in tasks:
        name = task_name(window, rep)
        for source, _ in _source_files(source_folded, window, rep):
            if not source.is_file():
                raise FileNotFoundError(f"pilot source is incomplete: {source}")
            if source.stat().st_size == 0:
                raise ValueError(f"pilot source file is empty: {source}")
            total_bytes += source.stat().st_size
        npz = source_folded / f"{name}.npz"
        _validate_baseline_npz(npz, pilot, window)
        key = f"fep/F64A/folded/{name}.npz"
        if key not in known_npz or sha256_file(npz) != known_npz[key]:
            raise ValueError(f"{npz}: checksum does not match the frozen forensic manifest")
        _validate_baseline_mdp(source_folded / name / "prod.mdp", window)
        reconcile_energy_history(
            source_folded / name / "dhdl.xvg", npz, cfg, pilot,
            float(pilot["baseline_total_time_ps"]),
        )
        stamp = (source_folded / name / "protocol.sha").read_text().strip()
        if stamp != pilot["source_protocol"]:
            raise ValueError(f"{name} source protocol is {stamp!r}")
    return {"tasks": len(tasks), "files": len(tasks) * 8, "bytes": total_bytes,
            "energy_links_verified": len(tasks)}


def stage(source_folded: str | Path, config: str | Path = DEFAULT_PILOT_CONFIG,
          first_light: bool = False) -> Path:
    """Copy and hash the immutable inputs needed by the continuation arm.

    Production staging requires all 60 checkpoints.  First-light staging copies only
    w0/r0 into its own throwaway namespace.  An existing stage is verified rather than
    silently replaced.
    """
    _, pilot = load_pilot(config)
    identity = _require_committed_code()
    source_folded = Path(source_folded).resolve()
    mode = "first_light" if first_light else "production"
    root = Path(pilot["_output_root"]) / mode
    staged = root / "source" / "F64A" / "folded"
    manifest_path = root / "stage_manifest.json"

    tasks = [(0, 0)] if first_light else [
        (w, r) for w in range(pilot["lambda_windows"])
        for r in range(pilot["replicates"])
    ]
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("pilot_id") != pilot["pilot_id"] or manifest.get("mode") != mode:
            raise ValueError(f"existing stage manifest belongs to another arm: {manifest_path}")
        if manifest.get("code_identity") != identity:
            raise ValueError(
                f"existing stage used different code identity: {manifest.get('code_identity')}"
            )
        if manifest.get("pilot_config_sha256") != sha256_file(pilot["_config_path"]):
            raise ValueError("pilot config changed after this stage was created")
        if manifest.get("base_config_sha256") != sha256_file(pilot["_base_config_path"]):
            raise ValueError("base pipeline config changed after this stage was created")
        for entry in manifest["files"]:
            path = root / entry["path"]
            if not path.is_file() or sha256_file(path) != entry["sha256"]:
                raise ValueError(f"existing stage failed checksum validation: {path}")
        print(f"Verified existing {mode} stage: {manifest_path}")
        return root
    if root.exists():
        raise FileExistsError(
            f"{root} exists without a stage manifest; refusing to reuse or overwrite it"
        )

    # Preflight the entire inventory before creating the destination. A missing late
    # checkpoint must not leave a plausible-looking partial stage behind.
    validate_source(source_folded, config, first_light)

    entries = []
    for window, rep in tasks:
        for source, relative in _source_files(source_folded, window, rep):
            destination = staged / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            _copy_atomic(source, destination)
            entries.append({
                "path": destination.relative_to(root).as_posix(),
                "size_bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
            })
        _validate_baseline_npz(staged / f"{task_name(window, rep)}.npz", pilot, window)
        _validate_baseline_mdp(staged / task_name(window, rep) / "prod.mdp", window)

    manifest = {
        "schema_version": 1,
        "pilot_id": pilot["pilot_id"],
        "mode": mode,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_path_recorded_for_audit": str(source_folded),
        "source_archive_was_modified": False,
        "pilot_config_sha256": sha256_file(pilot["_config_path"]),
        "base_config_sha256": sha256_file(pilot["_base_config_path"]),
        "code_identity": identity,
        "tasks": len(tasks),
        "files": sorted(entries, key=lambda row: row["path"]),
    }
    manifest_tmp = manifest_path.with_name(manifest_path.name + ".tmp")
    manifest_tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    manifest_tmp.replace(manifest_path)
    print(f"Staged {len(tasks)} checkpoint(s) into {root}")
    return root


def _load_stage_manifest(root: Path, pilot: dict, require_current_code: bool) -> dict:
    path = root / "stage_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(f"stage inputs before running: {path}")
    manifest = json.loads(path.read_text())
    if manifest.get("pilot_id") != pilot["pilot_id"]:
        raise ValueError(f"{path}: belongs to another pilot")
    if manifest.get("pilot_config_sha256") != sha256_file(pilot["_config_path"]):
        raise ValueError(f"{path}: pilot config changed after staging")
    if manifest.get("base_config_sha256") != sha256_file(pilot["_base_config_path"]):
        raise ValueError(f"{path}: base pipeline config changed after staging")
    identity = manifest.get("code_identity")
    if not isinstance(identity, dict) or not identity.get("git_commit"):
        raise ValueError(f"{path}: no code identity")
    if require_current_code:
        current = _require_committed_code()
        if current != identity:
            raise ValueError(
                f"{path}: staged under {identity}, but the current code is {current}; "
                "do not pull or edit pilot code while an arm is active"
            )
    return manifest


def _extension_protocol(pilot: dict, extension_ps: float) -> str:
    fields = {
        "pilot_id": pilot["pilot_id"],
        "continuation": "checkpoint_rng_and_coordinates",
        "source_protocol": pilot["source_protocol"],
        "extension_ps": extension_ps,
        "trajectory_retention": pilot["trajectory_retention"],
    }
    return hashlib.sha256(json.dumps(fields, sort_keys=True).encode()).hexdigest()[:16]


def _run_extension_mdrun(cfg: dict, run_dir: Path) -> None:
    """Run/resume the continuation, retrying only scheduler/GPU placement failures."""
    argv = [*mdrun_argv(cfg), "-s", "extension.tpr", "-deffnm", "prod",
            "-cpi", "prod.cpt", "-append", "-cpo", "prod.cpt",
            "-cpt", str(cfg["fep"]["checkpoint_interval_min"]),
            "-dhdl", "dhdl.xvg"]
    attempts = int(cfg["cluster"].get("gpu_retry_attempts", 3))
    delay = float(cfg["cluster"].get("gpu_retry_delay_s", 60))
    for attempt in range(1, attempts + 1):
        try:
            _run(argv, cwd=run_dir,
                 env={"OMP_NUM_THREADS": str(cfg["cluster"]["mdrun_ntomp"])})
            return
        except ToolError as exc:
            if not any(signature in str(exc) for signature in _GPU_BUSY_SIGNATURES):
                raise
            if attempt == attempts:
                raise GpuUnavailableError(
                    f"no usable GPU after {attempts} attempts; reschedule this task\n{exc}"
                ) from exc
            time.sleep(delay)


def _validate_output_npz(path: Path, pilot: dict, extension_ps: float,
                         stage_manifest: dict) -> None:
    """Fail unless one pilot output has the complete registered identity and shape."""
    expected_extension = int(round(extension_ps / float(pilot["energy_stride_ps"])))
    expected_shape = (
        int(pilot["lambda_windows"]),
        int(pilot["baseline_samples_per_window"]) + expected_extension,
    )
    expected_protocol = _extension_protocol(pilot, extension_ps)
    expected_identity = stage_manifest["code_identity"]
    with np.load(path) as data:
        required = {
            "u_kn_window", "lambda_index", "n_states", "provenance", "protocol",
            "pilot_id", "source_protocol", "source_npz_sha256", "baseline_samples",
            "extension_samples", "git_commit", "pilot_code_sha256",
            "pilot_config_sha256",
            "base_config_sha256",
        }
        missing = sorted(required - set(data.files))
        if missing:
            raise ValueError(f"{path}: missing output metadata {missing}")
        problems = []
        u = np.asarray(data["u_kn_window"])
        if u.shape != expected_shape:
            problems.append(f"shape={u.shape}, expected={expected_shape}")
        if not np.isfinite(u).all():
            problems.append("non-finite energies")
        expected_window = int(path.stem.split("_", 1)[0][1:])
        if int(data["lambda_index"]) != expected_window or int(data["n_states"]) != 20:
            problems.append("lambda metadata mismatch")
        expected_values = {
            "provenance": "gromacs_pmx",
            "protocol": expected_protocol,
            "pilot_id": pilot["pilot_id"],
            "source_protocol": pilot["source_protocol"],
            "git_commit": expected_identity["git_commit"],
            "pilot_code_sha256": expected_identity["pilot_code_sha256"],
            "pilot_config_sha256": stage_manifest["pilot_config_sha256"],
            "base_config_sha256": stage_manifest["base_config_sha256"],
        }
        source_key = f"source/F64A/folded/{path.name}"
        staged_hashes = {row["path"]: row["sha256"] for row in stage_manifest["files"]}
        if source_key not in staged_hashes:
            problems.append(f"stage manifest has no {source_key}")
        else:
            expected_values["source_npz_sha256"] = staged_hashes[source_key]
        for key, expected in expected_values.items():
            if str(data[key]) != str(expected):
                problems.append(f"{key}={data[key]!s}, expected={expected}")
        if int(data["baseline_samples"]) != pilot["baseline_samples_per_window"]:
            problems.append(f"baseline_samples={data['baseline_samples']}")
        if int(data["extension_samples"]) != expected_extension:
            problems.append(f"extension_samples={data['extension_samples']}")
    if problems:
        raise ValueError(f"{path}: invalid pilot output: {'; '.join(problems)}")


def run_task(window: int, rep: int, config: str | Path = DEFAULT_PILOT_CONFIG,
             first_light: bool = False) -> Path:
    """Continue one archived window and emit a baseline+extension NPZ."""
    cfg, pilot = load_pilot(config)
    _validate_task(pilot, window, rep)
    mode = "first_light" if first_light else "production"
    extension_ps = (float(pilot["first_light_extension_ps"]) if first_light
                    else float(pilot["extension_ns"]) * 1000.0)
    root = Path(pilot["_output_root"]) / mode
    stage_manifest = _load_stage_manifest(root, pilot, require_current_code=True)
    staged = root / "source" / "F64A" / "folded"
    name = task_name(window, rep)
    source_npz = staged / f"{name}.npz"
    source_run = staged / name
    _validate_baseline_npz(source_npz, pilot, window)

    out_dir = root / "fep" / "F64A" / "folded"
    out = out_dir / f"{name}.npz"
    if out.exists():
        _validate_output_npz(out, pilot, extension_ps, stage_manifest)
        print(f"Verified existing output {out}")
        return out

    run_dir = out_dir / name
    run_dir.mkdir(parents=True, exist_ok=True)
    protocol = _extension_protocol(pilot, extension_ps)
    stamp = run_dir / "extension_protocol.sha"
    if stamp.exists() and stamp.read_text().strip() != protocol:
        raise ToolError(f"{run_dir}: continuation protocol changed; refusing resume")
    if not stamp.exists() and (run_dir / "prod.cpt").exists():
        raise ToolError(f"{run_dir}: unstamped continuation checkpoint; refusing resume")
    stamp.write_text(protocol + "\n")

    # Materialize a private, checksum-identical append set once. The archive remains
    # read-only; subsequent scheduler invocations resume these local copies.
    for filename in ("prod.cpt", "prod.log", "prod.edr", "dhdl.xvg"):
        destination = run_dir / filename
        if not destination.exists():
            _copy_atomic(source_run / filename, destination)

    tpr = run_dir / "extension.tpr"
    if not tpr.exists():
        # GROMACS infers the output type from the final suffix and appends ``.tpr``
        # when it does not recognise one.  Keep the temporary name ending in .tpr
        # and pass its absolute path so the file we validate is exactly the file
        # convert-tpr writes.
        temporary_tpr = run_dir / f"extension.building.{os.getpid()}.tpr"
        _run([gmx_command(cfg), "convert-tpr", "-s",
              str(source_run / "prod.tpr"), "-extend", str(extension_ps),
              "-o", str(temporary_tpr)], cwd=run_dir)
        if not temporary_tpr.is_file() or temporary_tpr.stat().st_size == 0:
            raise ToolError(f"convert-tpr did not produce a usable {temporary_tpr}")
        temporary_tpr.replace(tpr)
    _run_extension_mdrun(cfg, run_dir)

    ext_path = run_dir / "dhdl.xvg"
    expected_end_ps = float(pilot["baseline_total_time_ps"]) + extension_ps
    u_ext = reconcile_energy_history(
        ext_path, source_npz, cfg, pilot, expected_end_ps,
    )
    expected_ext = (int(round(extension_ps / float(pilot["energy_stride_ps"]))))
    if u_ext.shape != (20, expected_ext):
        raise ValueError(
            f"{ext_path}: retained extension shape {u_ext.shape}, expected (20, {expected_ext})"
        )
    with np.load(source_npz) as baseline:
        combined = np.concatenate([np.asarray(baseline["u_kn_window"]), u_ext], axis=1)
    expected_total = pilot["baseline_samples_per_window"] + expected_ext
    if combined.shape != (20, expected_total) or not np.isfinite(combined).all():
        raise ValueError(f"{name}: invalid combined energy matrix {combined.shape}")

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.stem + ".tmp.npz")
    np.savez(
        tmp,
        lambda_index=window,
        n_states=20,
        u_kn_window=combined,
        provenance="gromacs_pmx",
        protocol=protocol,
        pilot_id=pilot["pilot_id"],
        source_protocol=pilot["source_protocol"],
        source_npz_sha256=sha256_file(source_npz),
        baseline_samples=pilot["baseline_samples_per_window"],
        extension_samples=expected_ext,
        git_commit=stage_manifest["code_identity"]["git_commit"],
        pilot_code_sha256=stage_manifest["code_identity"]["pilot_code_sha256"],
        pilot_config_sha256=stage_manifest["pilot_config_sha256"],
        base_config_sha256=stage_manifest["base_config_sha256"],
    )
    tmp.replace(out)
    _validate_output_npz(out, pilot, extension_ps, stage_manifest)
    print(f"Wrote {out} with {combined.shape[1]} samples/window")
    return out


def check_first_light(config: str | Path = DEFAULT_PILOT_CONFIG) -> Path:
    """Fail unless the isolated real-GPU continuation produced the expected result."""
    _, pilot = load_pilot(config)
    root = Path(pilot["_output_root"]) / "first_light"
    manifest = _load_stage_manifest(root, pilot, require_current_code=True)
    path = (root / "fep" / "F64A" /
            "folded" / "w0_r0.npz")
    if not path.is_file():
        raise FileNotFoundError(f"required first-light result is missing: {path}")
    _validate_output_npz(
        path, pilot, float(pilot["first_light_extension_ps"]), manifest,
    )
    print(f"First light passed: {path} (20, 3011), all finite")
    return path


def validate_production_outputs(
    config: str | Path = DEFAULT_PILOT_CONFIG,
) -> tuple[Path, dict]:
    """Validate identity, protocol, provenance, and exact shape of all 60 outputs."""
    _, pilot = load_pilot(config)
    root = Path(pilot["_output_root"]) / "production"
    manifest = _load_stage_manifest(root, pilot, require_current_code=False)
    fep_dir = root / "fep"
    extension_ps = float(pilot["extension_ns"]) * 1000.0
    for window in range(pilot["lambda_windows"]):
        for rep in range(pilot["replicates"]):
            path = fep_dir / "F64A" / "folded" / f"w{window}_r{rep}.npz"
            if not path.is_file():
                raise FileNotFoundError(f"missing production output: {path}")
            _validate_output_npz(path, pilot, extension_ps, manifest)
    print("Validated 60 production NPZs: one protocol/code identity, shape (20, 9001)")
    return fep_dir, manifest


def _solve_prefix(cfg: dict, pilot: dict, fep_dir: Path, rep: int,
                  stop: int, label: str, start: int = 0) -> dict:
    windows = []
    for window in range(pilot["lambda_windows"]):
        with np.load(fep_dir / "F64A" / "folded" / f"w{window}_r{rep}.npz") as data:
            raw = np.asarray(data["u_kn_window"])[:, start:stop]
        windows.append(decorrelate_window(raw, window))
    n_k = np.asarray([u.shape[1] for u in windows])
    solved = solve_leg_mbar(np.concatenate(windows, axis=1), n_k)
    kt = _KB_KCAL * float(cfg["fep"]["temperature_K"])
    local = adjacent_discrepancies(windows, kt)
    return {
        "replicate": rep,
        "checkpoint": label,
        "raw_column_start": start,
        "raw_column_stop": stop,
        "n_raw": int((stop - start) * pilot["lambda_windows"]),
        "n_independent": int(n_k.sum()),
        "dg_kcal": solved["dg_kT"] * kt,
        "mbar_error_kcal": solved["ddg_kT"] * kt,
        "net_hysteresis_kcal": float(abs(local.sum())),
        "sum_absolute_hysteresis_kcal": float(abs(local).sum()),
        "max_absolute_hysteresis_kcal": float(abs(local).max()),
        "min_adjacent_overlap": solved["min_adjacent"],
        "solver_notes": solved["solver_notes"],
    }


def analyze(config: str | Path = DEFAULT_PILOT_CONFIG) -> Path:
    """Analyze the 3, 6, and 9 ns folded-leg checkpoints without re-gating."""
    cfg, pilot = load_pilot(config)
    root = Path(pilot["_output_root"]) / "production"
    fep_dir, stage_manifest = validate_production_outputs(config)
    checkpoints = [(3001, "3_ns"), (6001, "6_ns"), (9001, "9_ns")]
    rows = [
        _solve_prefix(cfg, pilot, fep_dir, rep, stop, label)
        for stop, label in checkpoints
        for rep in range(pilot["replicates"])
    ]
    late = [
        _solve_prefix(cfg, pilot, fep_dir, rep, 9001, "final_3_ns", start=6001)
        for rep in range(pilot["replicates"])
    ]
    late_by_rep = {row["replicate"]: row for row in late}
    for row in rows:
        if row["checkpoint"] in {"3_ns", "6_ns"}:
            row["movement_vs_final_3_ns_kcal"] = abs(
                row["dg_kcal"] - late_by_rep[row["replicate"]]["dg_kcal"]
            )
            row["comparison_is_disjoint"] = True
        else:
            row["movement_vs_final_3_ns_kcal"] = None
            row["comparison_is_disjoint"] = False

    summaries = []
    for _, label in checkpoints:
        group = [row for row in rows if row["checkpoint"] == label]
        values = np.asarray([row["dg_kcal"] for row in group])
        propagated = np.sqrt(sum(row["mbar_error_kcal"] ** 2 for row in group)) / len(group)
        sem = float(np.std(values, ddof=1) / np.sqrt(len(values)))
        summaries.append({
            "checkpoint": label,
            "mean_dg_kcal": float(values.mean()),
            "uncertainty_kcal": float(max(propagated, sem)),
            "replicate_range_kcal": float(np.ptp(values)),
            "replicate_sem_kcal": sem,
        })
    result = {
        "schema_version": 1,
        "pilot_id": pilot["pilot_id"],
        "protocol": _extension_protocol(pilot, float(pilot["extension_ns"]) * 1000.0),
        "code_identity": stage_manifest["code_identity"],
        "pilot_config_sha256": stage_manifest["pilot_config_sha256"],
        "base_config_sha256": stage_manifest["base_config_sha256"],
        "scope": "F64A folded leg only; no experimental-accuracy or gate claim",
        "records": rows,
        "late_block_records": late,
        "summaries": summaries,
    }
    out = root / "analysis" / "f64a_extension.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(f"refusing to overwrite immutable pilot analysis {out}")
    tmp = out.with_name(out.name + ".tmp")
    tmp.write_text(json.dumps(result, indent=2) + "\n")
    tmp.replace(out)
    print(f"Wrote {out}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(DEFAULT_PILOT_CONFIG))
    sub = parser.add_subparsers(dest="command", required=True)
    p_stage = sub.add_parser("stage")
    p_stage.add_argument("--source-folded", required=True, type=Path)
    p_stage.add_argument("--first-light", action="store_true")
    p_validate = sub.add_parser("validate-source")
    p_validate.add_argument("--source-folded", required=True, type=Path)
    p_validate.add_argument("--first-light", action="store_true")
    sub.add_parser("verify-tools")
    sub.add_parser("check-first-light")
    sub.add_parser("validate-outputs")
    p_run = sub.add_parser("run")
    p_run.add_argument("--window", required=True, type=int)
    p_run.add_argument("--rep", required=True, type=int)
    p_run.add_argument("--first-light", action="store_true")
    sub.add_parser("analyze")
    args = parser.parse_args()

    try:
        if args.command == "stage":
            stage(args.source_folded, args.config, args.first_light)
        elif args.command == "validate-source":
            result = validate_source(args.source_folded, args.config, args.first_light)
            print(f"Validated {result['tasks']} tasks / {result['files']} files / "
                  f"{result['bytes']} bytes; reconciled "
                  f"{result['energy_links_verified']} XVG→NPZ energy histories")
        elif args.command == "verify-tools":
            cfg, _ = load_pilot(args.config)
            for key, value in verify_tools(cfg).items():
                print(f"{key}: {value}")
        elif args.command == "check-first-light":
            check_first_light(args.config)
        elif args.command == "validate-outputs":
            validate_production_outputs(args.config)
        elif args.command == "run":
            run_task(args.window, args.rep, args.config, args.first_light)
        else:
            analyze(args.config)
    except GpuUnavailableError as exc:
        print(f"RESCHEDULE: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(99) from exc


if __name__ == "__main__":
    main()
