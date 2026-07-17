"""Stage 1 prep helpers -- the parts that don't require OpenMM.

The heavy PDBFixer/Modeller path (:func:`prepare_variant`) is exercised as an
end-to-end run in the milestone, not here; these guard the pure logic.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.prep.build import (
    chains_for_oligomer,
    load_variant_record,
    mutation_string,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS_CSV = ROOT / "data" / "variants.csv"
STRUCTURE_3ECU = ROOT / "data" / "structures" / "3ECU.pdb"


def test_mutation_string():
    assert mutation_string("A", 4, "V") == "ALA-4-VAL"
    assert mutation_string("G", "93", "A") == "GLY-93-ALA"


def test_mutation_string_rejects_bad_code():
    with pytest.raises(KeyError):
        mutation_string("A", 4, "X")  # X is not a standard residue


def test_chains_monomer_needs_no_structure():
    assert chains_for_oligomer("does_not_exist.pdb", "monomer") == ["A"]


def test_chains_unknown_oligomer():
    with pytest.raises(ValueError, match="unknown oligomer"):
        chains_for_oligomer("x.pdb", "trimer")


@pytest.mark.skipif(not STRUCTURE_3ECU.exists(), reason="3ECU.pdb not fetched")
def test_chains_dimer_detects_partner():
    assert chains_for_oligomer(STRUCTURE_3ECU, "dimer") == ["A", "B"]


@pytest.mark.skipif(not VARIANTS_CSV.exists(), reason="variants.csv not built")
def test_load_variant_record():
    rec = load_variant_record(VARIANTS_CSV, "A4V")
    assert rec["wt_aa"] == "A" and rec["mut_aa"] == "V" and rec["oligomer"] == "monomer"
    with pytest.raises(KeyError):
        load_variant_record(VARIANTS_CSV, "Z999Z")
