"""FEP unfolded-leg reference builder (Stage 1).

Pure geometry helper is tested directly; the OpenMM/PDBFixer capping path is guarded on
the presence of both the toolkit and 3ECU (matching tests/test_prep.py convention) and
asserts the tripeptide parametrizes under amber14.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.prep.tripeptide import (
    _SKIP_RESNAMES,
    _unitvec,
    build_capped_tripeptide,
    rename_caps_for_gromacs,
)

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


class _StubAtom:
    def __init__(self, name):
        self.name = name


class _StubResidue:
    def __init__(self, name, atom_names):
        self.name = name
        self._atoms = [_StubAtom(n) for n in atom_names]

    def atoms(self):
        return iter(self._atoms)


class _StubTopology:
    def __init__(self, residues):
        self._residues = residues

    def residues(self):
        return iter(self._residues)


def _names(residue):
    return [a.name for a in residue.atoms()]


def test_rename_caps_maps_openmm_names_to_rtp():
    """pdb2gmx matches the amber .rtp: caps need HH31-33, and NME's methyl C is CH3."""
    ace = _StubResidue("ACE", ["CH3", "H1", "H2", "H3", "C", "O"])
    nme = _StubResidue("NME", ["N", "H", "C", "H1", "H2", "H3"])
    assert rename_caps_for_gromacs(_StubTopology([ace, nme])) == 7  # 3 in ACE, 4 in NME

    # ACE: methyl hydrogens renamed; the CARBONYL C and O must survive untouched.
    assert _names(ace) == ["CH3", "HH31", "HH32", "HH33", "C", "O"]
    # NME: PDBFixer renamed our CH3 to C -- put it back, or pdb2gmx aborts.
    assert _names(nme) == ["N", "H", "CH3", "HH31", "HH32", "HH33"]


def test_rename_caps_leaves_standard_residues_alone():
    """A backbone C in a normal residue must never be rewritten to CH3."""
    ala = _StubResidue("ALA", ["N", "H", "CA", "HA", "C", "O", "CB"])
    assert rename_caps_for_gromacs(_StubTopology([ala])) == 0
    assert _names(ala) == ["N", "H", "CA", "HA", "C", "O", "CB"]


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


def _angle_deg(p, q, r):
    v1, v2 = p - q, r - q
    return float(np.degrees(np.arccos(np.clip(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))))


def test_place_atom_reproduces_requested_internal_coordinates():
    """The NeRF fallback must build the geometry it is asked for, in the LOCAL frame.

    The construction it replaced put ACE's CH3, C and the residue's N and CA on one
    straight line (all angles 180 deg) and placed the carbonyl O along a fixed lab axis.
    """
    from src.prep.tripeptide import _place_atom

    a = np.array([0.0, 0.0, 0.0])
    b = np.array([0.15, 0.0, 0.0])
    c = np.array([0.22, 0.12, 0.0])
    d = _place_atom(a, b, c, 0.133, 121.7, 180.0)
    assert np.linalg.norm(d - c) == pytest.approx(0.133, abs=1e-6)
    assert _angle_deg(b, c, d) == pytest.approx(121.7, abs=1e-3)

    # equivariant under rotation of the reference frame -- no lab-axis dependence
    th = 0.7
    rot = np.array([[np.cos(th), -np.sin(th), 0.0],
                    [np.sin(th), np.cos(th), 0.0],
                    [0.0, 0.0, 1.0]])
    d_rot = _place_atom(a @ rot.T, b @ rot.T, c @ rot.T, 0.133, 121.7, 180.0)
    assert np.allclose(d_rot, d @ rot.T, atol=1e-9)


@pytest.mark.skipif(not (HAVE_OMM and STRUCTURE_3ECU.exists()),
                    reason="needs OpenMM/PDBFixer and 3ECU")
@pytest.mark.parametrize("center", [4, 93, 113])
def test_ace_cap_geometry_is_physical(center):
    """Cap heavy atoms must form a real peptide unit, not a collinear stub.

    ACE is chemically the preceding residue's CA-C(=O), so when that residue exists we
    copy its native coordinates. Previously every one of these angles was 180 degrees.
    """
    from openmm import unit

    top, positions = build_capped_tripeptide(str(STRUCTURE_3ECU), center, "A", 7.0)
    pos = np.array(positions.value_in_unit(unit.nanometer))
    at = {(a.residue.index, a.name): a.index for a in top.atoms()}
    residues = sorted({(a.residue.index, a.residue.name) for a in top.atoms()})
    ace, first = residues[0][0], residues[1][0]

    ch3, c, o = (pos[at[(ace, n)]] for n in ("CH3", "C", "O"))
    n_first, ca_first = pos[at[(first, "N")]], pos[at[(first, "CA")]]

    assert np.linalg.norm(c - n_first) == pytest.approx(0.133, abs=0.015)
    assert 105.0 < _angle_deg(ch3, c, n_first) < 130.0     # ideal ~116.6, was 180
    assert 105.0 < _angle_deg(c, n_first, ca_first) < 135.0  # ideal ~121.7, was 180
    assert 110.0 < _angle_deg(o, c, n_first) < 135.0       # ideal ~122.9, was arbitrary


@pytest.mark.skipif(not (HAVE_OMM and STRUCTURE_3ECU.exists()),
                    reason="needs OpenMM/PDBFixer and 3ECU")
def test_tripeptide_is_invariant_under_rotation_of_the_input(tmp_path):
    """The same protein in a different lab orientation must give the same molecule.

    The old cap construction referenced a fixed [0,0,1], so the carbonyl O moved relative
    to the peptide when the input PDB was rotated -- the unfolded leg was not reproducible.
    """
    from openmm import unit
    from openmm.app import PDBFile

    src = PDBFile(str(STRUCTURE_3ECU))
    pos = np.array(src.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    th = 0.9
    rot = np.array([[np.cos(th), -np.sin(th), 0.0],
                    [np.sin(th), np.cos(th), 0.0],
                    [0.0, 0.0, 1.0]])
    rotated = tmp_path / "rotated.pdb"
    with open(rotated, "w") as fh:
        PDBFile.writeFile(src.topology, unit.Quantity(pos @ rot.T, unit.nanometer), fh)

    def cap_angles(path):
        top, p = build_capped_tripeptide(str(path), 93, "A", 7.0)
        q = np.array(p.value_in_unit(unit.nanometer))
        at = {(a.residue.index, a.name): a.index for a in top.atoms()}
        res = sorted({(a.residue.index, a.residue.name) for a in top.atoms()})
        ace, first = res[0][0], res[1][0]
        return (_angle_deg(q[at[(ace, "CH3")]], q[at[(ace, "C")]], q[at[(first, "N")]]),
                _angle_deg(q[at[(ace, "O")]], q[at[(ace, "C")]], q[at[(first, "N")]]))

    for a, b in zip(cap_angles(STRUCTURE_3ECU), cap_angles(rotated)):
        assert a == pytest.approx(b, abs=0.5)
