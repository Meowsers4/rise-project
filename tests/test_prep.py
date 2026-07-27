"""Stage 1 prep helpers, plus the two invariants that must never silently break:
the disulfide really being reduced, and mature numbering really matching the PDB.

The heavy PDBFixer/Modeller path (:func:`prepare_variant`) is exercised as an
end-to-end run in the milestone, not here; these guard the pure logic. The two
invariant tests DO load the structure (guarded on OpenMM + 3ECU) because both
failures are invisible in pure logic -- they only show up in a built topology.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.prep.build import (
    assert_disulfide_reduced,
    assert_wt_residue,
    chains_for_oligomer,
    load_variant_record,
    mutation_string,
    strip_disulfide_bonds,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS_CSV = ROOT / "data" / "variants.csv"
STRUCTURE_3ECU = ROOT / "data" / "structures" / "3ECU.pdb"

try:
    import openmm  # noqa: F401
    import pdbfixer  # noqa: F401
    HAVE_OMM = True
except Exception:
    HAVE_OMM = False

needs_omm = pytest.mark.skipif(
    not (HAVE_OMM and STRUCTURE_3ECU.exists()),
    reason="OpenMM/PDBFixer or 3ECU.pdb unavailable",
)


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


# --------------------------------------------------------------------------- #
# Invariant 1: apo/reduced means REDUCED (README §2.1, CLAUDE.md #1)            #
# --------------------------------------------------------------------------- #
@needs_omm
def test_3ecu_ships_an_intact_disulfide():
    """Guards the premise of the fix: the crystal structure IS oxidized, so protonating
    it without breaking the bond leaves Cys57/Cys146 without HG."""
    from openmm.app import PDBFile
    top = PDBFile(str(STRUCTURE_3ECU)).topology
    ss = [(b[0].residue.id, b[1].residue.id) for b in top.bonds()
          if b[0].name == "SG" and b[1].name == "SG"]
    assert ("57", "146") in [tuple(sorted(p, key=int)) for p in ss]
    with pytest.raises(ValueError, match="NOT reduced"):
        assert_disulfide_reduced(top)


@needs_omm
def test_strip_disulfide_bonds_then_protonate_gives_free_thiols():
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=str(STRUCTURE_3ECU))
    fixer.removeChains([i for i, c in enumerate(fixer.topology.chains()) if c.id != "A"])
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    removed = strip_disulfide_bonds(fixer.topology)
    assert removed, "expected to remove the Cys57-Cys146 bond"
    fixer.addMissingHydrogens(7.0)
    assert_disulfide_reduced(fixer.topology)   # raises if any Cys lacks HG


# --------------------------------------------------------------------------- #
# Invariant 2: mature position addresses the residue we think it does (rule #3) #
# --------------------------------------------------------------------------- #
@needs_omm
def test_assert_wt_residue_accepts_mature_numbering():
    from openmm.app import PDBFile
    top = PDBFile(str(STRUCTURE_3ECU)).topology
    assert_wt_residue(top, ["A"], 4, "A")      # 3ECU resSeq 4 is ALA
    assert_wt_residue(top, ["A"], 93, "G")     # ... 93 is GLY


@needs_omm
def test_assert_wt_residue_catches_offset_numbering():
    """A precursor-numbered structure would shift every residue by mature_offset; the
    check must catch that instead of mutating whatever happens to sit there."""
    from openmm.app import PDBFile
    top = PDBFile(str(STRUCTURE_3ECU)).topology
    with pytest.raises(ValueError, match="is GLY, panel says ALA"):
        assert_wt_residue(top, ["A"], 93, "A")
    with pytest.raises(ValueError, match="no residue numbered"):
        assert_wt_residue(top, ["A"], 9999, "A")


@pytest.mark.skipif(not VARIANTS_CSV.exists(), reason="variants.csv not built")
@needs_omm
def test_every_panel_variant_matches_the_structure():
    """The whole panel must address real wild-type residues in the starting PDB."""
    import csv

    from openmm.app import PDBFile
    top = PDBFile(str(STRUCTURE_3ECU)).topology
    with open(VARIANTS_CSV) as f:
        for row in csv.DictReader(f):
            assert_wt_residue(top, ["A"], int(row["mature_pos"]), row["wt_aa"])
