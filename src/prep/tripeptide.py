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

from src.prep.build import load_config, load_variant_record

ROOT = Path(__file__).resolve().parents[2]

# Ideal bond lengths (nm) for placing cap heavy atoms; refined by minimisation downstream.
_D_PEPTIDE_CN = 0.133   # C(=O)--N peptide bond
_D_C_CH3 = 0.150        # carbonyl C--methyl C
_D_C_O = 0.123          # C=O
_D_N_CH3 = 0.147        # amide N--methyl C
_SKIP_RESNAMES = frozenset({"HOH", "WAT", "CU", "ZN", "NA", "CL"})


def _unitvec(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n else v


def build_capped_tripeptide(wt_pdb: str | Path, center: int, chain_id: str | None,
                            ph: float):
    """Return (topology, positions) for ACE-[i-1][i][i+1]-NME with hydrogens added.

    ``center`` is the residue id (structure numbering) of the mutation site; its two
    sequence neighbours are taken from the same chain. ``chain_id`` selects the protein
    chain (``None`` -> the chain that contains ``center``; first match wins -- monomer scope).
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
                break
        if trio:
            break
    if trio is None:
        raise KeyError(f"residue id {center} not found (chain={chain_id!r}) in {wt_pdb}")

    def coord(old_atom):
        return allpos[atom_index[old_atom]]

    def named(res, name):
        return next(a for a in res.atoms() if a.name == name)

    # --- cap heavy-atom geometry from the trio's terminal backbone atoms ---
    n0, ca0 = coord(named(trio[0], "N")), coord(named(trio[0], "CA"))
    d1 = _unitvec(n0 - ca0)
    c_ace = n0 + _D_PEPTIDE_CN * d1
    ch3_ace = c_ace + _D_C_CH3 * d1
    perp = _unitvec(np.cross(d1, np.array([0.0, 0.0, 1.0])))
    o_ace = c_ace + _D_C_O * perp

    cL, caL = coord(named(trio[-1], "C")), coord(named(trio[-1], "CA"))
    d2 = _unitvec(cL - caL)
    n_nme = cL + _D_PEPTIDE_CN * d2
    ch3_nme = n_nme + _D_N_CH3 * d2

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


def build_tripeptide(cfg: dict, variant: str, wt_pdb: str | Path, out_path: str | Path,
                     validate: bool = True) -> Path:
    """Write the unsolvated capped tripeptide reference for ``variant`` to ``out_path``."""
    from openmm.app import ForceField, PDBFile

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    center = int(record["mature_pos"])
    ph = cfg["structure"]["ph"]
    top, positions = build_capped_tripeptide(wt_pdb, center, chain_id=None, ph=ph)

    if validate:  # fail loudly here rather than deep inside a GPU job
        scfg = cfg["structure"]
        ff = ForceField("amber14-all.xml", f"amber14/{scfg['water_model']}.xml")
        ff.createSystem(top)

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
