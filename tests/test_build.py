"""Panel-builder logic (Stage 0), without network or the real structure.

Checks field derivation, the D4 oligomer rules (controls pinned to measured state;
uncharacterized via the geometric cutoff), and the controls floor -- all on a
committed synthetic fixture and temp CSVs, so it runs anywhere.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.panel.build import (
    _check_controls_floor,
    _oligomer_for_control,
    build_control_rows,
    build_uncharacterized_rows,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PDB = ROOT / "tests" / "fixtures" / "mini_dimer.pdb"

CFG = {
    "structure": {
        "default_oligomer": "monomer",
        "interface_cutoff_A": 5.0,
    },
    "panel": {"require_controls": True, "min_positive_controls": 2},
}


def test_oligomer_pinned_to_measured_state():
    assert _oligomer_for_control("apo_monomer", "monomer") == "monomer"
    assert _oligomer_for_control("apo_reduced", "monomer") == "monomer"
    assert _oligomer_for_control("holo_dimer", "monomer") == "dimer"
    # blank/unknown -> config default, never a guessed dimer
    assert _oligomer_for_control("", "monomer") == "monomer"


def test_build_control_rows_passthrough(tmp_path):
    csv_path = tmp_path / "controls.csv"
    csv_path.write_text(
        "variant,exp_ddg,exp_ddg_err,measured_state,exp_source,bucket\n"
        "A4V,NEEDS_VERIFICATION,NEEDS_VERIFICATION,apo_monomer,Lindberg2005,positive_control\n"
    )
    (row,) = build_control_rows(csv_path, CFG)
    assert (row["mature_pos"], row["wt_aa"], row["mut_aa"]) == (4, "A", "V")
    assert row["bucket"] == "positive_control"
    assert row["oligomer"] == "monomer"           # pinned to measured apo_monomer
    assert row["exp_ddg"] == "NEEDS_VERIFICATION"  # never fabricated
    assert row["exp_source"] == "Lindberg2005"


def test_build_uncharacterized_uses_geometry(tmp_path):
    # In the fixture, residue 1 is at the interface and residue 2 is buried.
    csv_path = tmp_path / "unchar.csv"
    csv_path.write_text(
        "variant,bucket,clinvar_id\n"
        "A1V,vus,VCV1\n"   # mature pos 1 -> interface -> dimer
        "G2A,vus,VCV2\n"   # mature pos 2 -> buried   -> monomer
    )
    rows = build_uncharacterized_rows(csv_path, CFG, FIXTURE_PDB)
    by_variant = {r["variant"]: r for r in rows}
    assert by_variant["A1V"]["oligomer"] == "dimer"
    assert by_variant["G2A"]["oligomer"] == "monomer"
    assert by_variant["A1V"]["exp_ddg"] == ""  # uncharacterized -> no experimental value


def test_controls_floor_enforced():
    rows = [{"bucket": "positive_control"}]  # only 1, floor is 2
    with pytest.raises(ValueError, match="positive controls"):
        _check_controls_floor(rows, CFG)
