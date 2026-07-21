"""FEP unfolded-leg reference builder (Stage 1).

Pure geometry helper is tested directly; the OpenMM/PDBFixer capping path is guarded on
the presence of both the toolkit and 3ECU (matching tests/test_prep.py convention) and
asserts the tripeptide parametrizes under amber14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.prep.tripeptide import _SKIP_RESNAMES, _unitvec, build_capped_tripeptide

ROOT = Path(__file__).resolve().parents[1]
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


def test_unitvec_normalizes():
    assert np.isclose(np.linalg.norm(_unitvec(np.array([3.0, 0.0, 4.0]))), 1.0)


def test_unitvec_zero_is_safe():
    assert np.allclose(_unitvec(np.zeros(3)), 0.0)  # no division-by-zero


def _first_protein_resid(pdb_path):
    from openmm.app import PDBFile
    top = PDBFile(str(pdb_path)).topology
    chain = next(c for c in top.chains() if c.id == "A")
    prot = [r for r in chain.residues() if r.name not in _SKIP_RESNAMES]
    return int(prot[0].id)


@needs_omm
def test_capped_tripeptide_parametrizes():
    from openmm.app import ForceField
    top, _ = build_capped_tripeptide(STRUCTURE_3ECU, center=50, chain_id="A", ph=7.0)
    names = [r.name for r in top.residues()]
    assert names[0] == "ACE" and names[-1] == "NME"
    assert len(names) == 5                      # ACE + trio + NME
    assert names[2] == "PHE"                    # SOD1 residue 50 is the central residue
    # createSystem raises if any residue/atom fails to match an amber14 template
    ForceField("amber14-all.xml", "amber14/tip3p.xml").createSystem(top)


@needs_omm
def test_missing_center_raises():
    with pytest.raises(KeyError):
        build_capped_tripeptide(STRUCTURE_3ECU, center=99999, chain_id="A", ph=7.0)


@needs_omm
def test_terminus_center_raises():
    first = _first_protein_resid(STRUCTURE_3ECU)
    with pytest.raises(ValueError, match="terminus"):
        build_capped_tripeptide(STRUCTURE_3ECU, center=first, chain_id="A", ph=7.0)
