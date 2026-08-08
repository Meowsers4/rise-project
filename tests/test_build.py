"""Panel-builder logic (Stage 0), without network or the real structure.

Checks control ingestion (apo-monomer pinning, err default, dimer cross-check
passthrough), the Route B constraints (missense-only, mature range, bucket, dedupe),
oligomer-by-geometry, and the controls floor -- all on a committed synthetic fixture
and temp CSVs, so it runs anywhere.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.panel.build import (
    _check_controls_floor,
    _validate_gate_subset,
    build_control_rows,
    build_uncharacterized_rows,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PDB = ROOT / "tests" / "fixtures" / "mini_dimer.pdb"

CFG = {
    "project": {"n_residues_mature": 153},
    "structure": {"default_oligomer": "monomer", "interface_cutoff_A": 5.0},
    "panel": {"require_controls": True, "min_positive_controls": 2},
    "validation": {"gate_subset": ["A4V"]},
}

_CONTROLS_HEADER = (
    "variant,wt_aa,mature_pos,mut_aa,exp_ddg,exp_ddg_err,exp_ddg_dimer,class,oligomer,exp_source\n"
)


def test_control_ingest_pins_monomer_and_carries_dimer(tmp_path):
    csv_path = tmp_path / "controls.csv"
    csv_path.write_text(_CONTROLS_HEADER + "A4V,A,4,V,1.62,0.3,4.31,B,monomer,Kumar2017\n")
    (row,) = build_control_rows(csv_path, CFG)
    assert (row["mature_pos"], row["wt_aa"], row["mut_aa"]) == (4, "A", "V")
    assert row["bucket"] == "positive_control"
    assert row["oligomer"] == "monomer"          # pinned
    assert row["measured_state"] == "apo_monomer"
    assert row["exp_ddg"] == "1.62"              # positive = destabilizing, verbatim
    assert row["exp_ddg_dimer"] == "4.31"        # cross-check carried, not gated


def test_control_err_defaults_when_blank(tmp_path):
    csv_path = tmp_path / "controls.csv"
    csv_path.write_text(_CONTROLS_HEADER + "I18V,I,18,V,0.37,,,B,monomer,Kumar2017\n")
    (row,) = build_control_rows(csv_path, CFG)
    assert row["exp_ddg_err"] == "0.30"


def test_control_column_mismatch_errors(tmp_path):
    # variant says position 4 but the row tabulates 5 -> loud failure
    csv_path = tmp_path / "controls.csv"
    csv_path.write_text(_CONTROLS_HEADER + "A4V,A,5,V,1.0,0.3,,B,monomer,src\n")
    with pytest.raises(ValueError, match="disagrees"):
        build_control_rows(csv_path, CFG)


def test_uncharacterized_uses_geometry_and_dedupes(tmp_path):
    # Fixture: residue 1 at interface -> dimer; residue 2 buried -> monomer.
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text(
        "variant,bucket,clinvar_id\nA1V,vus,VCV1\nG2A,pathogenic_uncharacterized,VCV2\n"
    )
    rows = build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB, control_variants=set())
    by = {r["variant"]: r for r in rows}
    assert by["A1V"]["oligomer"] == "dimer"
    assert by["G2A"]["oligomer"] == "monomer"
    assert by["A1V"]["exp_ddg"] == ""  # no ground truth for payoff targets


@pytest.mark.parametrize("bad_row", ["G127X,vus,x", "E134*,vus,x", "A4,vus,x"])
def test_uncharacterized_rejects_non_missense(tmp_path, bad_row):
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text("variant,bucket,clinvar_id\n" + bad_row + "\n")
    with pytest.raises(ValueError):
        build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB, control_variants=set())


def test_uncharacterized_rejects_bad_bucket(tmp_path):
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text("variant,bucket,clinvar_id\nA1V,positive_control,x\n")
    with pytest.raises(ValueError, match="bucket"):
        build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB, control_variants=set())


def test_uncharacterized_dedupe_collision_errors(tmp_path):
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text("variant,bucket,clinvar_id\nA4V,vus,x\n")
    with pytest.raises(ValueError, match="collides"):
        build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB, control_variants={"A4V"})


def test_uncharacterized_out_of_range_errors(tmp_path):
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text("variant,bucket,clinvar_id\nA200V,vus,x\n")
    with pytest.raises(ValueError, match="out of range"):
        build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB, control_variants=set())


def test_gate_subset_must_be_controls():
    _validate_gate_subset({"A4V"}, CFG)  # ok
    with pytest.raises(ValueError, match="gate_subset"):
        _validate_gate_subset({"G93A"}, CFG)  # A4V missing from controls


def test_controls_floor_enforced():
    rows = [{"bucket": "positive_control"}]  # only 1, floor is 2
    with pytest.raises(ValueError, match="positive controls"):
        _check_controls_floor(rows, CFG)


def test_panel_rebuild_preserves_the_audit_columns():
    """rule panel regenerates variants.csv wholesale from OUTPUT_FIELDS.

    A column missing from the writer is silently deleted on the next rebuild. That is
    what happened to charge_change / wells2021 / axakova_class, which carry claim C1,
    the Wells-10 head-to-head, and claim C3 respectively.
    """
    from src.panel.build import OUTPUT_FIELDS, audit_columns

    for col in ("charge_change", "wells2021", "axakova_class"):
        assert col in OUTPUT_FIELDS, f"{col} would be dropped by a panel rebuild"

    # derived, not hand-maintained
    assert audit_columns("D90A", "D", "A")["charge_change"] == "true"
    assert audit_columns("A4V", "A", "V")["charge_change"] == "false"
    assert audit_columns("H46R", "H", "R")["charge_change"] == "true"   # His = neutral
    assert audit_columns("A4V", "A", "V")["wells2021"] == "true"
    assert audit_columns("G93V", "G", "V")["wells2021"] == "false"


def test_committed_panel_is_reproducible_from_its_sources():
    """data/variants.csv must equal what a rebuild produces.

    Otherwise a Snakemake rebuild silently changes the panel underneath the gate.
    """
    import tempfile

    from src.panel.build import build_panel, load_config

    cfg = load_config(ROOT / "config" / "pipeline.yaml")
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "variants.csv"
        build_panel(cfg, out)
        assert out.read_text() == (ROOT / "data" / "variants.csv").read_text(), \
            "data/variants.csv is stale relative to its build inputs"


def test_gate_subset_rejects_charge_changing_variants():
    """CLAUDE.md rule 2, enforced in code instead of by a comment.

    Net-charge PME artifacts do not cancel between the folded and unfolded legs and the
    pmx engine has no counterion co-alchemy, so a charge-changing variant in the gate
    would be judged on an uncancelled few-kcal/mol error.
    """
    from src.panel.build import _validate_gate_subset

    rows = [
        {"variant": "A4V", "charge_change": "false"},
        {"variant": "H43R", "charge_change": "true"},
    ]
    controls = {"A4V", "H43R"}

    ok = {"validation": {"gate_subset": ["A4V"]}}
    _validate_gate_subset(controls, ok, rows)          # neutral only: fine

    bad = {"validation": {"gate_subset": ["A4V", "H43R"]}}
    with pytest.raises(ValueError, match="charge-CHANGING"):
        _validate_gate_subset(controls, bad, rows)
