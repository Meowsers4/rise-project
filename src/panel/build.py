"""Stage 0 -- build and stratify the SOD1 variant panel.

Produces ``data/variants.csv`` (the panel source of truth) with the schema from
README §4:

    variant,mature_pos,wt_aa,mut_aa,bucket,exp_ddg,exp_source,clinvar_id,oligomer

Inputs:
- ``data/controls_experimental.csv`` -- curated positive controls with experimental
  ΔΔG (the M3 validation set). Values may be ``NEEDS_VERIFICATION`` until filled from
  literature; they are passed through verbatim so nothing is fabricated here.
- (optional) ``data/uncharacterized_variants.csv`` -- ClinVar/ALSoD/gnomAD-derived
  variants with columns ``variant,bucket,clinvar_id``. Until the external ingestion
  is implemented (see :func:`ingest_external_variants`), the panel is controls-only.

Oligomer assignment follows D4:
- **Controls** are pinned to the oligomeric state their ΔΔG was MEASURED in
  (``measured_state``), never the geometric rule -- this protects the gate.
- **Uncharacterized** variants use the geometric 5 Å interface rule
  (:func:`src.panel.interface.classify_oligomer`) against ``structure.starting_pdb``.

All parameters (offset, cutoff, default oligomer, starting PDB, controls floor) are
read from ``config/pipeline.yaml``; none are hardcoded here (CLAUDE.md rule 4).
"""

from __future__ import annotations

import argparse
import csv
import urllib.request
from pathlib import Path

import yaml

from src.panel.interface import classify_oligomer
from src.panel.numbering import parse_variant

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_FIELDS = [
    "variant",
    "mature_pos",
    "wt_aa",
    "mut_aa",
    "bucket",
    "exp_ddg",
    "exp_ddg_err",
    "exp_source",
    "measured_state",
    "clinvar_id",
    "oligomer",
]

# measured_state values that mean the control was characterised as a single subunit.
_MONOMER_STATES = {"apo_monomer", "apo_reduced", "monomer"}
_DIMER_STATES = {"apo_dimer", "holo_dimer", "dimer"}


def load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def ensure_structure(cfg: dict) -> Path:
    """Return the path to ``structure.starting_pdb``, downloading it if absent.

    Structures live under ``data/structures/`` (gitignored). Fetched from RCSB.
    """
    pdb_id = cfg["structure"]["starting_pdb"]
    dest = ROOT / "data" / "structures" / f"{pdb_id}.pdb"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, dest)  # noqa: S310 (trusted RCSB host)
    return dest


def _oligomer_for_control(measured_state: str, default: str) -> str:
    """Pin a control to the oligomeric state it was measured in (D4 caveat)."""
    state = measured_state.strip().lower()
    if state in _MONOMER_STATES:
        return "monomer"
    if state in _DIMER_STATES:
        return "dimer"
    # Unknown/blank measured state -> fall back to config default, do not guess a dimer.
    return default


def _base_row(variant: str) -> dict:
    wt, pos, mut = parse_variant(variant)
    return {"variant": variant, "mature_pos": pos, "wt_aa": wt, "mut_aa": mut}


def build_control_rows(controls_csv: str | Path, cfg: dict) -> list[dict]:
    """Rows for the curated experimental controls (oligomer pinned to measured state)."""
    default_oligomer = cfg["structure"]["default_oligomer"]
    rows: list[dict] = []
    with open(controls_csv) as f:
        for entry in csv.DictReader(f):
            variant = entry["variant"].strip()
            row = _base_row(variant)
            row.update(
                bucket=entry.get("bucket", "positive_control") or "positive_control",
                exp_ddg=entry.get("exp_ddg", ""),
                exp_ddg_err=entry.get("exp_ddg_err", ""),
                exp_source=entry.get("exp_source", ""),
                measured_state=entry.get("measured_state", ""),
                clinvar_id=entry.get("clinvar_id", ""),
                oligomer=_oligomer_for_control(
                    entry.get("measured_state", ""), default_oligomer
                ),
            )
            rows.append(row)
    return rows


def build_uncharacterized_rows(
    uncharacterized_csv: str | Path, cfg: dict, pdb_path: str | Path
) -> list[dict]:
    """Rows for uncharacterized variants; oligomer via the geometric 5 Å rule (D4)."""
    cutoff = cfg["structure"]["interface_cutoff_A"]
    default_oligomer = cfg["structure"]["default_oligomer"]
    rows: list[dict] = []
    with open(uncharacterized_csv) as f:
        for entry in csv.DictReader(f):
            variant = entry["variant"].strip()
            row = _base_row(variant)
            row.update(
                bucket=entry.get("bucket", "vus") or "vus",
                exp_ddg="",
                exp_ddg_err="",
                exp_source="",
                measured_state="",
                clinvar_id=entry.get("clinvar_id", ""),
                oligomer=classify_oligomer(
                    row["mature_pos"], str(pdb_path), cutoff, default=default_oligomer
                ),
            )
            rows.append(row)
    return rows


def ingest_external_variants(cfg: dict) -> Path:
    """Pull ClinVar + ALSoD, cross-ref gnomAD -> ``data/uncharacterized_variants.csv``.

    Not yet implemented -- this is the external-data stage. Rather than fabricate a
    panel, we raise so the caller knows to supply the CSV. See README §4 Stage 0.
    """
    raise NotImplementedError(
        "ClinVar/ALSoD/gnomAD ingestion is not implemented yet. Provide "
        "data/uncharacterized_variants.csv (columns: variant,bucket,clinvar_id) "
        "or run with controls only."
    )


def build_panel(cfg: dict, controls_csv: str | Path, out_path: str | Path) -> list[dict]:
    """Assemble the panel (controls + any uncharacterized) and write ``variants.csv``."""
    pdb_path = ensure_structure(cfg)
    rows = build_control_rows(controls_csv, cfg)

    uncharacterized_csv = ROOT / "data" / "uncharacterized_variants.csv"
    if uncharacterized_csv.exists():
        rows += build_uncharacterized_rows(uncharacterized_csv, cfg, pdb_path)

    _check_controls_floor(rows, cfg)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _check_controls_floor(rows: list[dict], cfg: dict) -> None:
    """Fail fast if the panel has fewer positive controls than the config floor."""
    if not cfg["panel"].get("require_controls", True):
        return
    n_controls = sum(1 for r in rows if r["bucket"] == "positive_control")
    floor = cfg["panel"].get("min_positive_controls", 1)
    if n_controls < floor:
        raise ValueError(
            f"Panel has {n_controls} positive controls; config requires >= {floor}. "
            "Add controls to data/controls_experimental.csv."
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the SOD1 variant panel (Stage 0).")
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--controls", default=str(ROOT / "data" / "controls_experimental.csv"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    rows = build_panel(cfg, args.controls, args.out)
    n_controls = sum(1 for r in rows if r["bucket"] == "positive_control")
    print(f"Wrote {len(rows)} variants ({n_controls} positive controls) -> {args.out}")


if __name__ == "__main__":
    main()
