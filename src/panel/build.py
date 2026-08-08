"""Stage 0 -- build and stratify the SOD1 variant panel.

Produces ``data/variants.csv`` (the panel source of truth) with the schema from
README §4:

    variant,mature_pos,wt_aa,mut_aa,bucket,exp_ddg,exp_source,clinvar_id,oligomer

Inputs (paths from ``config.panel``):
- ``panel.controls_csv`` (``data/sod1_controls_apo_monomer.csv``) -- verified apo-monomer
  positive controls (Kumar2017 Table S1). ``exp_ddg`` is the apo-monomer gate value
  (positive = destabilizing); ``exp_ddg_dimer`` is a carried cross-check only. Values
  pass through verbatim -- nothing is fabricated here.
- (optional) ``panel.uncharacterized_csv`` (``data/uncharacterized_variants.csv``) --
  Route B input with columns ``variant,bucket,clinvar_id``. Until the external
  ingestion is implemented (see :func:`ingest_external_variants`) this file is
  user-supplied; if it has no data rows the panel is controls-only.

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
    "exp_ddg",        # apo-monomer ddG; positive = destabilizing (the gate value)
    "exp_ddg_err",
    "exp_ddg_dimer",  # dimer cross-check only -- NOT used in the monomer gate
    "exp_source",
    "measured_state",
    "clinvar_id",
    "oligomer",
    # Added by the 2026-08-07 literature audit (README §4 Stage 0). These MUST stay in the
    # writer: rule panel regenerates variants.csv wholesale, so a column missing here is
    # silently deleted on the next rebuild, taking claim C1 (charge_change), the Wells-10
    # head-to-head (wells2021) and C3 (axakova_class) with it.
    "charge_change",
    "wells2021",
    "axakova_class",
]

# Default experimental error (kcal/mol) when a controls row leaves exp_ddg_err blank
# (Kumar2017 Table S1 stated methodological error).
_DEFAULT_DDG_ERR = "0.30"

# --- audit columns (README §4 Stage 0, added 2026-08-07) ------------------------------
# Wells NGM et al., Protein Science 2021;30(9):1804 -- the 10 charge-conserving variants
# that study covered. Our head-to-head comparison set (claim C2 / milestone M3.5).
_WELLS_2021 = frozenset(
    ("A4V", "A4S", "A4T", "L38V", "G41S", "G93A", "G93S", "I113T", "V148G", "V148I")
)
# Formal side-chain charge at pH 7. His is NEUTRAL here (pKa ~6.0), which matches what the
# pipeline actually builds: PDBFixer at structure.ph 7.0 and pdb2gmx both default every
# histidine to a neutral tautomer. Changing this changes which variants claim C1 covers.
_SIDECHAIN_CHARGE = {"D": -1, "E": -1, "K": +1, "R": +1}


def audit_columns(variant: str, wt_aa: str, mut_aa: str) -> dict[str, str]:
    """The three columns the literature audit requires, DERIVED not hand-maintained.

    Kept in one place so the controls and VUS paths cannot disagree, and so a panel
    rebuild reproduces them instead of dropping them.
    """
    dq = _SIDECHAIN_CHARGE.get(mut_aa, 0) - _SIDECHAIN_CHARGE.get(wt_aa, 0)
    return {
        "charge_change": "true" if dq else "false",
        "wells2021": "true" if variant in _WELLS_2021 else "false",
        # Populated from the Axakova 2025 DMS supplementary; blocks claim C3 until then.
        "axakova_class": "",
    }




# Buckets permitted for Route B uncharacterized variants.
_UNCHARACTERIZED_BUCKETS = {"pathogenic_uncharacterized", "vus"}


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


def _parse_and_check(variant: str, entry: dict) -> tuple[str, int, str]:
    """parse_variant + cross-check against any wt_aa/mature_pos/mut_aa in the row.

    Rejects non-missense names (parse_variant raises on ``G127X``/``E134*``/etc.) and
    catches a variant string that disagrees with its own tabulated columns.
    """
    wt, pos, mut = parse_variant(variant)
    for col, want, got in (
        ("wt_aa", entry.get("wt_aa"), wt),
        ("mut_aa", entry.get("mut_aa"), mut),
        ("mature_pos", entry.get("mature_pos"), pos),
    ):
        if want not in (None, ""):
            if str(want).strip() != str(got):
                raise ValueError(
                    f"{variant}: column {col}={want!r} disagrees with parsed {got!r}"
                )
    return wt, pos, mut


def build_control_rows(controls_csv: str | Path, cfg: dict) -> list[dict]:
    """Rows for the verified apo-monomer controls.

    Oligomer is pinned to ``monomer`` (measured state), overriding the 5 Å rule (D4).
    ``exp_ddg`` (apo-monomer, positive=destabilizing) is the gate value; ``exp_ddg_dimer``
    is carried as a cross-check only. Values pass through verbatim -- no fabrication.
    """
    rows: list[dict] = []
    with open(controls_csv) as f:
        for entry in csv.DictReader(f):
            variant = entry["variant"].strip()
            wt, pos, mut = _parse_and_check(variant, entry)
            file_oligomer = (entry.get("oligomer") or "monomer").strip().lower()
            if file_oligomer != "monomer":
                raise ValueError(
                    f"{variant}: controls are apo-monomer measurements but oligomer="
                    f"{file_oligomer!r}; refusing to mislabel the pinned state."
                )
            rows.append(
                {
                    "variant": variant,
                    "mature_pos": pos,
                    "wt_aa": wt,
                    "mut_aa": mut,
                    "bucket": "positive_control",
                    "exp_ddg": entry["exp_ddg"].strip(),
                    "exp_ddg_err": (entry.get("exp_ddg_err") or "").strip() or _DEFAULT_DDG_ERR,
                    "exp_ddg_dimer": (entry.get("exp_ddg_dimer") or "").strip(),
                    "exp_source": (entry.get("exp_source") or "").strip(),
                    "measured_state": "apo_monomer",
                    "clinvar_id": (entry.get("clinvar_id") or "").strip(),
                    "oligomer": "monomer",  # pinned; overrides geometry (D4)
                    **audit_columns(variant, wt, mut),
                }
            )
    return rows


def build_uncharacterized_rows(
    uncharacterized_csv: str | Path,
    cfg: dict,
    pdb_path: str | Path,
    control_variants: set[str],
) -> list[dict]:
    """Rows for Route B uncharacterized variants, with constraints enforced loudly.

    Missense-only (via parse_variant), mature numbering in range, bucket in the allowed
    set, and NO overlap with the controls. Oligomer comes from the geometric 5 Å rule.
    """
    cutoff = cfg["structure"]["interface_cutoff_A"]
    default_oligomer = cfg["structure"]["default_oligomer"]
    n_mature = cfg["project"]["n_residues_mature"]
    rows: list[dict] = []
    seen: set[str] = set()
    with open(uncharacterized_csv) as f:
        for entry in csv.DictReader(f):
            variant = entry["variant"].strip()
            bucket = (entry.get("bucket") or "").strip()
            if bucket not in _UNCHARACTERIZED_BUCKETS:
                raise ValueError(
                    f"{variant}: bucket {bucket!r} not in {sorted(_UNCHARACTERIZED_BUCKETS)}"
                )
            wt, pos, mut = _parse_and_check(variant, entry)  # rejects non-missense
            if not 1 <= pos <= n_mature:
                raise ValueError(
                    f"{variant}: mature position {pos} out of range 1..{n_mature}"
                )
            if variant in control_variants:
                raise ValueError(f"{variant}: collides with a control variant (dedupe)")
            if variant in seen:
                raise ValueError(f"{variant}: duplicated in {uncharacterized_csv}")
            seen.add(variant)
            rows.append(
                {
                    "variant": variant,
                    "mature_pos": pos,
                    "wt_aa": wt,
                    "mut_aa": mut,
                    "bucket": bucket,
                    "exp_ddg": "",
                    "exp_ddg_err": "",
                    "exp_ddg_dimer": "",
                    "exp_source": "",
                    "measured_state": "",
                    "clinvar_id": (entry.get("clinvar_id") or "").strip(),
                    "oligomer": classify_oligomer(
                        pos, str(pdb_path), cutoff, default=default_oligomer
                    ),
                    **audit_columns(variant, wt, mut),
                }
            )
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


def _validate_gate_subset(control_variants: set[str], cfg: dict,
                          rows: list[dict] | None = None) -> None:
    """Every M3 gate-subset variant must be a charge-neutral control.

    CLAUDE.md rule 2 makes charge-neutrality non-negotiable: a net-charge change carries a
    PME finite-size artifact that does not cancel between the folded and unfolded legs
    (different box sizes), and the pmx engine has no counterion co-alchemy. Until now that
    was guaranteed only by the hand-curated list in config and a comment above it, so
    adding H43R to lift `n` would have passed every check.
    """
    subset = cfg.get("validation", {}).get("gate_subset") or []
    missing = [v for v in subset if v not in control_variants]
    if missing:
        raise ValueError(
            f"validation.gate_subset variants not in controls: {missing}"
        )
    if rows is not None:
        charged = sorted(r["variant"] for r in rows
                         if r["variant"] in subset and r.get("charge_change") == "true")
        if charged:
            raise ValueError(
                f"validation.gate_subset contains charge-CHANGING variants {charged}. "
                "The gate is charge-neutral only (CLAUDE.md rule 2): net-charge PME "
                "artifacts do not cancel between legs and the engine has no counterion "
                "co-alchemy. These belong in the C1 sub-gate, not this one."
            )


def build_panel(cfg: dict, out_path: str | Path) -> list[dict]:
    """Assemble the panel (controls + any uncharacterized) and write ``variants.csv``.

    Input paths come from ``config.panel`` (``controls_csv``, ``uncharacterized_csv``).
    The uncharacterized file is optional; if it has no data rows the panel is controls
    only. Uncharacterized variants are deduped against the controls.
    """
    pdb_path = ensure_structure(cfg)

    controls_csv = ROOT / cfg["panel"]["controls_csv"]
    rows = build_control_rows(controls_csv, cfg)
    control_variants = {r["variant"] for r in rows}
    _validate_gate_subset(control_variants, cfg, rows)

    uncharacterized_csv = ROOT / cfg["panel"]["uncharacterized_csv"]
    if uncharacterized_csv.exists():
        rows += build_uncharacterized_rows(
            uncharacterized_csv, cfg, pdb_path, control_variants
        )

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
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    rows = build_panel(cfg, args.out)
    n_controls = sum(1 for r in rows if r["bucket"] == "positive_control")
    n_unchar = len(rows) - n_controls
    olig = {}
    for r in rows:
        olig[r["oligomer"]] = olig.get(r["oligomer"], 0) + 1
    print(f"Wrote {len(rows)} variants -> {args.out}")
    print(f"  controls: {n_controls}  uncharacterized: {n_unchar}")
    print(f"  oligomer: {dict(sorted(olig.items()))}")


if __name__ == "__main__":
    main()
