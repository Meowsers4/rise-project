"""Stage 1 (FEP unfolded-leg reference) -- build the capped tripeptide for a variant.

The unfolded state of the folding thermodynamic cycle is modelled as a solvated capped
tripeptide carrying the mutation site with its NATIVE SOD1 neighbours (decision locked
2026-07-21): ``ACE-[i-1][i][i+1]-NME``. The three residue identities and starting
coordinates are sliced from the WT prepared structure (single source of truth -- no
separately stored sequence), capped with ACE/NME, and written UNSOLVATED; Perses solvates
and alchemically mutates the central residue in :mod:`src.fep.window`.

Method (validated against amber14 ``createSystem``): build the chain in physical order
ACE -> trio -> NME (Modeller/PDBFixer use residue order for terminus detection), place the
cap heavy atoms by ideal peptide-bond geometry, then let PDBFixer add all hydrogens. Cap
geometry is a starting guess; Perses minimises before sampling.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from src.prep.build import _THREE, load_config, load_variant_record

ROOT = Path(__file__).resolve().parents[2]

# Ideal bond lengths (nm) for placing cap heavy atoms; refined by minimisation downstream.
_D_PEPTIDE_CN = 0.133   # C(=O)--N peptide bond
_D_C_CH3 = 0.150        # carbonyl C--methyl C
_D_C_O = 0.123          # C=O
_D_N_CH3 = 0.147        # amide N--methyl C
_SKIP_RESNAMES = frozenset({"HOH", "WAT", "CU", "ZN", "NA", "CL"})

# PDBFixer names the caps by OpenMM convention: methyl hydrogens H1/H2/H3, and it RENAMES
# the NME methyl carbon we create as "CH3" to "C". GROMACS' amber .rtp entries use
# HH31/HH32/HH33 and keep the NME methyl carbon as CH3, so pdb2gmx aborts with
# "Atom C in residue NME 5 was not found in rtp entry NME with 6 atoms" (verified against
# amber99sb-star-ildn-mut, 2026-08-01). ACE's carbonyl C and O already match and are NOT
# renamed -- only NME's C is the methyl. Applied to the WRITTEN file only; the in-memory
# topology keeps OpenMM names so ForceField.createSystem() still validates it.
_CAP_ATOM_RENAMES: dict[str, dict[str, str]] = {
    "ACE": {"H1": "HH31", "H2": "HH32", "H3": "HH33"},
    "NME": {"C": "CH3", "H1": "HH31", "H2": "HH32", "H3": "HH33"},
}


# Ideal backbone angles (degrees) for the fallback construction. Only reached when the
# mutation site is two residues from a chain terminus, i.e. never for the current panel.
_A_C_N_CA = 121.7       # ACE C -- res0 N -- res0 CA
_A_CH3_C_N = 116.6      # ACE CH3 -- ACE C -- res0 N
_A_O_C_N = 122.9        # ACE O -- ACE C -- res0 N
_A_CA_C_N = 116.2       # resN CA -- resN C -- NME N
_A_C_N_CH3 = 121.7      # resN C -- NME N -- NME CH3
_OMEGA_TRANS = 180.0    # trans peptide bond


def _unitvec(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n else v


def _place_atom(a: np.ndarray, b: np.ndarray, c: np.ndarray,
                bond: float, angle_deg: float, dihedral_deg: float) -> np.ndarray:
    """Position a fourth atom D from internal coordinates against A-B-C (NeRF).

    ``bond`` is |C-D| (nm), ``angle_deg`` is B-C-D, ``dihedral_deg`` is A-B-C-D. Built in
    the local frame of A-B-C, so the result rotates with the input rather than depending
    on a lab axis -- the defect this replaced placed cap atoms along a fixed [0,0,1].
    """
    ang, dih = np.radians(angle_deg), np.radians(dihedral_deg)
    bc = _unitvec(c - b)
    n = _unitvec(np.cross(b - a, bc))
    local = np.array([
        -bond * np.cos(ang),
        bond * np.sin(ang) * np.cos(dih),
        bond * np.sin(ang) * np.sin(dih),
    ])
    basis = np.array([bc, np.cross(n, bc), n]).T
    return c + basis @ local


def build_capped_tripeptide(wt_pdb: str | Path, center: int, chain_id: str | None,
                            ph: float, expect_resname: str | None = None):
    """Return (topology, positions) for ACE-[i-1][i][i+1]-NME with hydrogens added.

    ``center`` is the residue id (structure numbering) of the mutation site; its two
    sequence neighbours are taken from the same chain. ``chain_id`` selects the protein
    chain (``None`` -> the chain that contains ``center``; first match wins -- monomer scope).
    ``expect_resname`` is the three-letter wild-type residue the panel says sits at
    ``center``; when given, a mismatch raises rather than silently building the unfolded
    reference for the wrong residue.
    """
    from openmm import Vec3, unit
    from openmm.app import PDBFile, Topology
    from openmm.app.element import carbon, nitrogen, oxygen
    from pdbfixer import PDBFixer

    pdb = PDBFile(str(wt_pdb))
    top = pdb.topology
    allpos = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    atom_index = {a: i for i, a in enumerate(top.atoms())}

    # locate the mutation-site residue and its two flanking residues (same chain)
    def protein_residues(chain):
        return [r for r in chain.residues() if r.name not in _SKIP_RESNAMES]

    trio = None
    for chain in top.chains():
        if chain_id is not None and chain.id != chain_id:
            continue
        residues = protein_residues(chain)
        for i, r in enumerate(residues):
            if int(r.id) == center:
                if i == 0 or i == len(residues) - 1:
                    raise ValueError(
                        f"residue {center} is at a chain terminus in {wt_pdb}; cannot build "
                        "a tripeptide with both native neighbours."
                    )
                trio = residues[i - 1:i + 2]
                # The cap heavy atoms ARE backbone atoms of the residues flanking the
                # trio, when those exist. Carrying them lets us copy native coordinates
                # instead of constructing a guess -- see _ace_cap_coords.
                flank_prev = residues[i - 2] if i - 2 >= 0 else None
                flank_next = residues[i + 2] if i + 2 < len(residues) else None
                # Guard the neighbours we slice: `residues` is list-adjacent after
                # filtering, not necessarily sequence-adjacent. Across a chain gap the
                # "neighbour" is a non-native residue with no bond to the trio.
                ids = [int(r.id) for r in trio]
                if ids != [center - 1, center, center + 1]:
                    raise ValueError(
                        f"residues flanking {center} in {wt_pdb} are numbered {ids}, not "
                        f"{[center - 1, center, center + 1]} -- there is a chain break, so "
                        "these are not the native neighbours the unfolded reference needs."
                    )
                break
        if trio:
            break
    if trio is None:
        raise KeyError(f"residue id {center} not found (chain={chain_id!r}) in {wt_pdb}")
    if expect_resname is not None and trio[1].name != expect_resname:
        raise ValueError(
            f"residue {center} in {wt_pdb} is {trio[1].name}, panel says {expect_resname}. "
            "Refusing to build the unfolded reference for the wrong residue."
        )

    def coord(old_atom):
        return allpos[atom_index[old_atom]]

    def named(res, name):
        return next(a for a in res.atoms() if a.name == name)

    # --- cap heavy-atom geometry ---------------------------------------------------
    # An ACE cap is chemically the preceding residue's CA-C(=O); an NME cap is the
    # following residue's N-CA. When those residues are present in the parent structure
    # -- which they are for every panel variant except one two residues from a terminus
    # -- copy their coordinates. That is exact native geometry and, unlike a constructed
    # guess, is invariant under rigid-body rotation of the input PDB.
    n0, ca0 = coord(named(trio[0], "N")), coord(named(trio[0], "CA"))
    cL, caL = coord(named(trio[-1], "C")), coord(named(trio[-1], "CA"))

    if flank_prev is not None:
        ch3_ace = coord(named(flank_prev, "CA"))
        c_ace = coord(named(flank_prev, "C"))
        o_ace = coord(named(flank_prev, "O"))
    else:
        c0 = coord(named(trio[0], "C"))
        c_ace = _place_atom(c0, ca0, n0, _D_PEPTIDE_CN, _A_C_N_CA, _OMEGA_TRANS)
        ch3_ace = _place_atom(ca0, n0, c_ace, _D_C_CH3, _A_CH3_C_N, _OMEGA_TRANS)
        o_ace = _place_atom(ca0, n0, c_ace, _D_C_O, _A_O_C_N, 0.0)

    if flank_next is not None:
        n_nme = coord(named(flank_next, "N"))
        ch3_nme = coord(named(flank_next, "CA"))
    else:
        nL = coord(named(trio[-1], "N"))
        n_nme = _place_atom(nL, caL, cL, _D_PEPTIDE_CN, _A_CA_C_N, _OMEGA_TRANS)
        ch3_nme = _place_atom(caL, cL, n_nme, _D_N_CH3, _A_C_N_CH3, _OMEGA_TRANS)

    # --- assemble ACE -> trio -> NME in physical order ---
    new = Topology()
    ch = new.addChain("A")
    pos: list[np.ndarray] = []

    ace = new.addResidue("ACE", ch, id="0")
    a_ch3 = new.addAtom("CH3", carbon, ace)
    a_c = new.addAtom("C", carbon, ace)
    a_o = new.addAtom("O", oxygen, ace)
    pos += [ch3_ace, c_ace, o_ace]
    new.addBond(a_ch3, a_c)
    new.addBond(a_c, a_o)

    gmap: dict = {}
    first_n = last_c = None
    for r in trio:
        nr = new.addResidue(r.name, ch, id=r.id)
        for a in r.atoms():
            gmap[a] = new.addAtom(a.name, a.element, nr)
            pos.append(coord(a))
            if r is trio[0] and a.name == "N":
                first_n = gmap[a]
            if r is trio[-1] and a.name == "C":
                last_c = gmap[a]
    for b in top.bonds():  # preserve intra- and inter-residue (backbone) connectivity
        if b[0] in gmap and b[1] in gmap:
            new.addBond(gmap[b[0]], gmap[b[1]])
    new.addBond(a_c, first_n)  # ACE.C -- res0.N

    nme = new.addResidue("NME", ch, id="999")
    n_n = new.addAtom("N", nitrogen, nme)
    n_ch3 = new.addAtom("CH3", carbon, nme)
    pos += [n_nme, ch3_nme]
    new.addBond(n_n, n_ch3)
    new.addBond(last_c, n_n)  # resN.C -- NME.N

    positions = unit.Quantity(np.array(pos), unit.nanometer)

    # --- add hydrogens via PDBFixer (handles ACE/NME caps; Modeller mis-detects termini) ---
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        tmp = Path(fh.name)
        PDBFile.writeFile(new, positions, fh)
    try:
        fixer = PDBFixer(filename=str(tmp))
        fixer.findMissingResidues()
        fixer.missingResidues = {}          # do NOT extend the capped termini
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)       # reduced Cys get HG; caps get methyl H's
    finally:
        tmp.unlink(missing_ok=True)
    return fixer.topology, fixer.positions


def rename_caps_for_gromacs(topology) -> int:
    """Rename ACE/NME cap atoms in place, OpenMM convention -> GROMACS amber ``.rtp``.

    Returns the number of atoms renamed (7 for a well-formed capped tripeptide: 3 methyl
    hydrogens in ACE, plus the methyl carbon and 3 hydrogens in NME). Call this
    only AFTER any OpenMM ``createSystem`` validation: ``amber14-all.xml`` matches the
    OpenMM names, and renaming first makes the template lookup fail. See
    :data:`_CAP_ATOM_RENAMES` for why the names differ.
    """
    renamed = 0
    for residue in topology.residues():
        mapping = _CAP_ATOM_RENAMES.get(residue.name)
        if not mapping:
            continue
        for atom in residue.atoms():
            new_name = mapping.get(atom.name)
            if new_name is not None:
                atom.name = new_name
                renamed += 1
    return renamed


def build_tripeptide(cfg: dict, variant: str, wt_pdb: str | Path, out_path: str | Path,
                     validate: bool = True) -> Path:
    """Write the unsolvated capped tripeptide reference for ``variant`` to ``out_path``."""
    from openmm.app import ForceField, PDBFile

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    center = int(record["mature_pos"])
    ph = cfg["structure"]["ph"]
    top, positions = build_capped_tripeptide(wt_pdb, center, chain_id=None, ph=ph,
                                             expect_resname=_THREE[record["wt_aa"]])

    if validate:  # fail loudly here rather than deep inside a GPU job
        scfg = cfg["structure"]
        ff = ForceField("amber14-all.xml", f"amber14/{scfg['water_model']}.xml")
        ff.createSystem(top)

    # Stage 3's pmx engine feeds this file to pdb2gmx, which matches atom names against
    # the amber .rtp. Rename after validation, so the file on disk is GROMACS-readable.
    rename_caps_for_gromacs(top)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        PDBFile.writeFile(top, positions, f)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the capped tripeptide FEP unfolded-leg reference (Stage 1).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--wt-pdb", required=True,
                        help="WT apo (unsolvated) protein from `build.py --wt`.")
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = build_tripeptide(cfg, args.variant, args.wt_pdb, args.out)
    print(f"Built capped tripeptide for {args.variant} -> {out}")


if __name__ == "__main__":
    main()
