"""Dimer-interface classification for SOD1 variants (D4).

A variant is simulated as a ``dimer`` if the mutated residue sits at the SOD1
homodimer interface, else as a ``monomer`` (``structure.default_oligomer``). "At the
interface" is defined geometrically -- any heavy atom of the residue within
``structure.interface_cutoff_A`` of the partner subunit -- so the residue set is
DERIVED, never a hardcoded list (D4 / CLAUDE.md rule 3).

The reference structure (``structure.starting_pdb``, 3ECU) uses mature numbering
(1-153), so a structure ``resSeq`` equals the mature variant position directly. The
apo crystal has two dimers in the asymmetric unit (chains A/B/C/D); the partner of a
reference chain is detected by contact count, not by assuming a fixed pairing.

Dependency-light on purpose: a numpy-only PDB reader, so the classification logic and
its test run without mdtraj/MDAnalysis installed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Generous radius used ONLY to identify which chain is the biological dimer partner
# (robustness heuristic, independent of the physics interface cutoff from config).
# The partner is unambiguous well beyond any sane interface cutoff, so this is not a
# tunable pipeline parameter and stays a documented internal constant.
_PARTNER_DETECT_CUTOFF_A = 8.0


@dataclass(frozen=True)
class _Atom:
    chain: str
    res_seq: int
    element: str
    xyz: tuple[float, float, float]


def _parse_pdb_heavy_atoms(pdb_path: str) -> list[_Atom]:
    """Read heavy (non-hydrogen) ATOM records from a PDB file.

    HETATM records (waters, ions) are ignored, so an apo structure yields only
    protein heavy atoms. Column offsets follow the PDB v3.3 fixed-width spec.
    """
    atoms: list[_Atom] = []
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            element = line[76:78].strip() or line[12:16].strip()[:1]
            if element == "H":
                continue
            atoms.append(
                _Atom(
                    chain=line[21],
                    res_seq=int(line[22:26]),
                    element=element,
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                )
            )
    if not atoms:
        raise ValueError(f"No ATOM records parsed from {pdb_path!r}")
    return atoms


def _chain_coords(atoms: list[_Atom], chain: str) -> np.ndarray:
    return np.array([a.xyz for a in atoms if a.chain == chain], dtype=float)


def _min_dist_per_residue(
    atoms: list[_Atom], ref_chain: str, partner_coords: np.ndarray
) -> dict[int, float]:
    """Per-residue minimum heavy-atom distance from ``ref_chain`` to ``partner_coords``."""
    out: dict[int, float] = {}
    ref = [a for a in atoms if a.chain == ref_chain]
    if partner_coords.size == 0:
        return out
    for atom in ref:
        d = float(np.min(np.linalg.norm(partner_coords - np.array(atom.xyz), axis=1)))
        prev = out.get(atom.res_seq)
        if prev is None or d < prev:
            out[atom.res_seq] = d
    return out


def detect_partner_chain(
    pdb_path: str, ref_chain: str, cutoff_A: float = _PARTNER_DETECT_CUTOFF_A
) -> str:
    """Return the chain forming the true dimer interface with ``ref_chain``.

    The partner is the other chain with the most heavy-atom contacts within
    ``cutoff_A`` -- robust to multiple dimers in the asymmetric unit. ``cutoff_A``
    defaults to a generous detection radius, decoupled from the (possibly tighter)
    physics interface cutoff, so partner detection never fails just because the
    interface cutoff is small.
    """
    atoms = _parse_pdb_heavy_atoms(pdb_path)
    chains = sorted({a.chain for a in atoms})
    if ref_chain not in chains:
        raise ValueError(f"Reference chain {ref_chain!r} not in {chains}")
    ref_coords = _chain_coords(atoms, ref_chain)
    best_chain, best_contacts = None, -1
    for chain in chains:
        if chain == ref_chain:
            continue
        other = _chain_coords(atoms, chain)
        # count ref atoms with any partner atom within cutoff
        contacts = int(
            np.sum(
                [
                    np.any(np.linalg.norm(other - p, axis=1) <= cutoff_A)
                    for p in ref_coords
                ]
            )
        )
        if contacts > best_contacts:
            best_chain, best_contacts = chain, contacts
    if best_chain is None or best_contacts <= 0:
        raise ValueError(
            f"No inter-chain contacts within {cutoff_A} A for chain {ref_chain!r}"
        )
    return best_chain


def interface_residues(
    pdb_path: str,
    cutoff_A: float,
    ref_chain: str = "A",
    partner_chain: str | None = None,
) -> set[int]:
    """Mature-numbered residues of ``ref_chain`` within ``cutoff_A`` of the partner.

    If ``partner_chain`` is None it is auto-detected via :func:`detect_partner_chain`.
    """
    if partner_chain is None:
        partner_chain = detect_partner_chain(pdb_path, ref_chain)
    atoms = _parse_pdb_heavy_atoms(pdb_path)
    partner_coords = _chain_coords(atoms, partner_chain)
    per_res = _min_dist_per_residue(atoms, ref_chain, partner_coords)
    return {res for res, dist in per_res.items() if dist <= cutoff_A}


def classify_oligomer(
    mature_pos: int,
    pdb_path: str,
    cutoff_A: float,
    default: str = "monomer",
    ref_chain: str = "A",
) -> str:
    """Return ``"dimer"`` if ``mature_pos`` is an interface residue, else ``default``.

    ``default`` should be ``structure.default_oligomer`` from config.
    """
    at_interface = interface_residues(pdb_path, cutoff_A, ref_chain=ref_chain)
    return "dimer" if mature_pos in at_interface else default
