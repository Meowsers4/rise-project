"""Create a deterministic checksum manifest for the forensic analysis archives.

The manifest records the raw NPZ windows without copying them into git, plus the JSON,
topology, and repository inputs consumed by the CPU-only reanalysis and figure builders.
Absolute machine-specific root paths are deliberately omitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of ``path`` using bounded memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _entry(root_id: str, root: Path, path: Path) -> dict:
    return {
        "root": root_id,
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def aggregate_digest(entries: list[dict]) -> str:
    """Return a stable digest over sorted file identity, size, and checksum fields."""
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda row: (row["root"], row["path"])):
        line = (
            f"{entry['root']}\0{entry['path']}\0{entry['size_bytes']}\0"
            f"{entry['sha256']}\n"
        )
        digest.update(line.encode())
    return digest.hexdigest()


def _selected_paths(attempt1_archive: Path, ss_archive: Path,
                    repository: Path) -> dict[str, list[Path]]:
    attempt1 = list(attempt1_archive.glob("fep/*/*/w*_r*.npz"))
    attempt1.extend(attempt1_archive.glob("fep/*/ddg.json"))
    attempt1.extend(attempt1_archive.glob("convergence/*.json"))
    attempt1.extend(attempt1_archive.glob("archive/A4V_3ns_*.json"))
    attempt1.extend(attempt1_archive.glob("reanalysis_2026-09-12/convergence/*.json"))
    attempt1.extend(attempt1_archive.glob("fep/G93A/folded/system_r*/hybrid.top"))

    ss = list(ss_archive.glob("*/w*_r*.npz"))
    ss.extend(ss_archive.glob("*.json"))
    ss.extend(ss_archive.glob("MANIFEST.md"))
    ss.extend(ss_archive.glob("folded/system_r*/hybrid.top"))

    repo = [
        repository / "config" / "pipeline.yaml",
        repository / "data" / "reference_sensitivity_scenarios.csv",
        repository / "src" / "analysis" / "reference_sensitivity.py",
        repository / "src" / "analysis" / "methods_figures.py",
        repository / "src" / "analysis" / "archive_manifest.py",
    ]
    selected = {
        "main_archive": sorted(set(attempt1)),
        "ss_archive": sorted(set(ss)),
        "repository": sorted(set(repo)),
    }
    for root_id, paths in selected.items():
        missing = [path for path in paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(f"missing {root_id} inputs: {missing}")
    return selected


def build_manifest(attempt1_archive: str | Path, ss_archive: str | Path,
                   repository: str | Path = ROOT, workers: int = 4) -> dict:
    """Hash the complete raw-window inventories and direct analysis inputs."""
    attempt1_archive = Path(attempt1_archive).resolve()
    ss_archive = Path(ss_archive).resolve()
    repository = Path(repository).resolve()
    roots = {
        "main_archive": attempt1_archive,
        "ss_archive": ss_archive,
        "repository": repository,
    }
    selected = _selected_paths(attempt1_archive, ss_archive, repository)

    cfg = load_config(repository / "config" / "pipeline.yaml")
    n_windows = int(cfg["fep"]["lambda_windows"])
    n_replicates = int(cfg["fep"]["replicates"])
    n_legs = 2
    n_variants = len(cfg["validation"]["gate_subset"])
    expected_attempt1_npz = n_variants * n_legs * n_windows * n_replicates
    expected_ss_npz = n_legs * n_windows * n_replicates
    actual_main_npz = sum(path.suffix == ".npz" for path in selected["main_archive"])
    actual_a4v_npz = sum(
        path.suffix == ".npz" and path.relative_to(attempt1_archive).parts[:2] == ("fep", "A4V")
        for path in selected["main_archive"]
    )
    actual_rederivable_attempt1_npz = actual_main_npz - actual_a4v_npz
    actual_ss_npz = sum(path.suffix == ".npz" for path in selected["ss_archive"])
    if actual_main_npz != expected_attempt1_npz:
        raise ValueError(
            f"main archive raw inventory is {actual_main_npz}, expected {expected_attempt1_npz}"
        )
    expected_one_variant = n_legs * n_windows * n_replicates
    if actual_a4v_npz != expected_one_variant:
        raise ValueError(f"A4V raw inventory is {actual_a4v_npz}, expected {expected_one_variant}")
    if actual_ss_npz != expected_ss_npz:
        raise ValueError(f"SS raw inventory is {actual_ss_npz}, expected {expected_ss_npz}")

    jobs = [
        (root_id, roots[root_id], path)
        for root_id, paths in selected.items()
        for path in paths
    ]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        entries = list(pool.map(lambda args: _entry(*args), jobs))
    entries.sort(key=lambda row: (row["root"], row["path"]))
    return {
        "schema_version": 1,
        "root_names": {root_id: root.name for root_id, root in roots.items()},
        "inventory": {
            "main_archive_raw_npz": actual_main_npz,
            "gate_attempt1_rederivable_raw_npz": actual_rederivable_attempt1_npz,
            "a4v_attempt2_raw_npz": actual_a4v_npz,
            "ss_raw_npz": actual_ss_npz,
            "files_total": len(entries),
            "bytes_total": sum(row["size_bytes"] for row in entries),
        },
        "aggregate_sha256": aggregate_digest(entries),
        "files": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt1-archive", required=True, type=Path)
    parser.add_argument("--ss-archive", required=True, type=Path)
    parser.add_argument("--repository", default=ROOT, type=Path)
    parser.add_argument("--output", default=ROOT / "data/forensic_archive_manifest.json",
                        type=Path)
    parser.add_argument("--workers", default=4, type=int)
    args = parser.parse_args()
    manifest = build_manifest(
        args.attempt1_archive,
        args.ss_archive,
        args.repository,
        args.workers,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"Wrote {args.output}: {manifest['inventory']['files_total']} files, "
        f"{manifest['inventory']['bytes_total']} bytes, "
        f"aggregate {manifest['aggregate_sha256']}"
    )


if __name__ == "__main__":
    main()
