"""Stage 1 -- prepare an apo, disulfide-reduced SOD1 system for one variant.

A parameterized ``variant_id -> prepared_system`` step (README §4 Stage 1), run once
per variant with no manual steps. Pipeline (PDBFixer + OpenMM ``Modeller``):

1. Load ``structure.starting_pdb`` (3ECU) and keep only the chains for this variant's
   oligomeric state (monomer -> ref chain; dimer -> ref + geometric partner).
2. Apply the point mutation in every kept chain (a homodimer carries it in both).
3. Rebuild 3ECU's disordered metal-binding loops (D5) via PDBFixer missing-residue/
   atom addition -- internal gaps only, not terminal tails.
4. Strip all heteroatoms (metals + crystal waters) -> apo (§2.1).
5. Break the Cys57-Cys146 disulfide, then protonate at ``structure.ph`` so both
   cysteines get HG (free thiols) -> disulfide REDUCED, and VERIFY it.
6. Solvate (``water_model``) with ``solvent_padding_nm`` padding and neutralize to
   ``ion_conc_M``.

Step 5 is not automatic. 3ECU ships ``SSBOND CYS A 57 - CYS A 146`` at 2.04 A; OpenMM
infers that bond on load and ``addMissingHydrogens`` then (correctly) refuses to add HG
to a bonded cysteine -- so simply protonating an oxidized crystal structure yields the
OXIDIZED form, silently, which is the opposite of the disease-relevant species this
project is built around. :func:`strip_disulfide_bonds` removes the bond first and
:func:`assert_disulfide_reduced` fails the build if the result is not actually reduced.
The freed thiols start ~2 A apart, i.e. clashing; energy minimisation before production
must relax them (``fep.minimize``).

Everything is read from ``config/pipeline.yaml``; nothing is hardcoded (CLAUDE.md #4).
Apo-first and reduced-disulfide are enforced here, not optional (§2.1, CLAUDE.md #1).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml

from src.panel.build import ensure_structure
from src.panel.interface import detect_partner_chain

ROOT = Path(__file__).resolve().parents[2]

# one-letter -> three-letter amino acid codes (PDBFixer applyMutations wants 3-letter).
_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def load_config(config_path: str | Path) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_variant_record(variants_csv: str | Path, variant: str) -> dict:
    """Return the panel row for ``variant`` (variants.csv is the source of truth)."""
    with open(variants_csv) as f:
        for row in csv.DictReader(f):
            if row["variant"] == variant:
                return row
    raise KeyError(f"{variant!r} not found in {variants_csv}")


def mutation_string(wt_aa: str, mature_pos: int | str, mut_aa: str) -> str:
    """PDBFixer mutation spec, e.g. ('A', 4, 'V') -> 'ALA-4-VAL'."""
    return f"{_THREE[wt_aa]}-{int(mature_pos)}-{_THREE[mut_aa]}"


def chains_for_oligomer(pdb_path: str | Path, oligomer: str, ref_chain: str = "A") -> list[str]:
    """Chain ids to keep: [ref] for monomer, [ref, partner] for dimer (geometry-derived)."""
    if oligomer == "monomer":
        return [ref_chain]
    if oligomer == "dimer":
        return [ref_chain, detect_partner_chain(str(pdb_path), ref_chain)]
    raise ValueError(f"unknown oligomer {oligomer!r} (expected monomer|dimer)")


def assert_wt_residue(topology, chain_ids: list[str], mature_pos: int, wt_aa: str) -> None:
    """Verify the residue at ``mature_pos`` really is ``wt_aa`` in every kept chain.

    The pipeline addresses residues by MATURE position and relies on the starting PDB
    using mature numbering (3ECU does: resSeq 1-153). Nothing else checks that, so a
    structure numbered any other way -- e.g. ``structure.fallback_pdb`` -- would mutate a
    silently wrong residue and every ΔΔG downstream would be for a different variant.

    Args:
        topology: OpenMM ``Topology`` to inspect.
        chain_ids: Chain ids that must each carry the residue.
        mature_pos: Mature-numbered residue position from the panel.
        wt_aa: Expected one-letter wild-type residue from the panel.

    Raises:
        ValueError: If the residue is absent or is not the expected wild type.
    """
    want = _THREE[wt_aa]
    for cid in chain_ids:
        chain = next((c for c in topology.chains() if c.id == cid), None)
        if chain is None:
            raise ValueError(f"chain {cid!r} not found in the starting structure")
        res = next((r for r in chain.residues() if r.id == str(mature_pos)), None)
        if res is None:
            raise ValueError(
                f"chain {cid}: no residue numbered {mature_pos}. The starting structure "
                "is not in mature numbering (project.mature_offset applies to the "
                "precursor, not to residue ids in the PDB)."
            )
        if res.name != want:
            raise ValueError(
                f"chain {cid}: residue {mature_pos} is {res.name}, panel says {want} "
                f"({wt_aa}). Refusing to mutate the wrong residue -- check that "
                "structure.starting_pdb uses mature (153-residue) numbering."
            )


def strip_disulfide_bonds(topology) -> list[tuple[str, str]]:
    """Remove every SG-SG bond from ``topology`` in place; return the pairs removed.

    Must run AFTER ``addMissingAtoms`` (which rebuilds the topology) and BEFORE
    ``addMissingHydrogens``: OpenMM decides CYS-vs-CYX from these bonds, so the bond has
    to be gone before protonation for the thiol hydrogens to be added.
    """
    removed = [
        (b[0].residue.id, b[1].residue.id)
        for b in topology.bonds()
        if b[0].name == "SG" and b[1].name == "SG"
    ]
    if removed:
        # Topology exposes no public bond removal; _bonds is the list bonds() iterates.
        topology._bonds = [
            b for b in topology._bonds
            if not (b[0].name == "SG" and b[1].name == "SG")
        ]
    return removed


def assert_disulfide_reduced(topology) -> None:
    """Fail unless every cysteine is a free thiol (has HG) and no SG-SG bond survives.

    Raises:
        ValueError: If the structure is oxidized -- v1 is reduced-only (CLAUDE.md #1).
    """
    surviving = [
        (b[0].residue.id, b[1].residue.id)
        for b in topology.bonds()
        if b[0].name == "SG" and b[1].name == "SG"
    ]
    unprotonated = [
        f"{r.name}{r.id}"
        for r in topology.residues()
        if any(a.name == "SG" for a in r.atoms()) and not any(a.name == "HG" for a in r.atoms())
    ]
    if surviving or unprotonated:
        raise ValueError(
            "disulfide is NOT reduced: "
            f"surviving SG-SG bonds={surviving}, cysteines without HG={unprotonated}. "
            "v1 simulates the apo, disulfide-REDUCED form (README §2.1, CLAUDE.md #1)."
        )


def _trim_terminal_missing(fixer) -> None:
    """Drop terminal missing residues so we rebuild internal loops only, not new tails."""
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]


def prepare_variant(cfg: dict, variant: str, out_path: str | Path,
                    apply_mutation: bool = True, solvate: bool = True) -> Path:
    """Build the apo, reduced system for ``variant`` and write it to PDB.

    Defaults (``apply_mutation=True, solvate=True``) produce the solvated MUTANT used
    by the FoldX/Rosetta prescreen (Stage 2). For the Perses FEP folded leg pass
    ``apply_mutation=False, solvate=False``: Perses needs the **WT** apo/reduced protein
    UNSOLVATED and builds the WT->mut hybrid + solvent itself (see
    :func:`prepare_wt_apo` and ``src/fep/window.py``).
    """
    # Imports are local so the pure helpers above (and their tests) don't need OpenMM.
    from openmm import unit
    from openmm.app import ForceField, Modeller, PDBFile
    from pdbfixer import PDBFixer

    scfg = cfg["structure"]
    if scfg["form"] != "apo" or scfg["disulfide"] != "reduced":
        raise ValueError(
            "v1 is apo + disulfide-reduced only (§2.1). "
            f"Got form={scfg['form']!r} disulfide={scfg['disulfide']!r} -- STOP and ask."
        )

    starting_pdb = ensure_structure(cfg)
    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    keep = chains_for_oligomer(starting_pdb, record["oligomer"])
    mut = mutation_string(record["wt_aa"], record["mature_pos"], record["mut_aa"])

    fixer = PDBFixer(filename=str(starting_pdb))

    # keep only the chains for this oligomeric state
    all_chain_ids = [c.id for c in fixer.topology.chains()]
    remove = [i for i, cid in enumerate(all_chain_ids) if cid not in keep]
    fixer.removeChains(remove)

    # the panel addresses residues by mature position; prove that matches this structure
    # BEFORE mutating anything (CLAUDE.md #3).
    assert_wt_residue(fixer.topology, keep, int(record["mature_pos"]), record["wt_aa"])

    # mutate every kept chain (homodimer carries the mutation in both subunits).
    # Skipped for the FEP folded leg -- Perses introduces the mutation alchemically.
    if apply_mutation:
        for cid in keep:
            fixer.applyMutations([mut], cid)

    # rebuild disordered internal loops (D5); or none if disabled
    fixer.findMissingResidues()
    if cfg["structure"].get("rebuild_missing_loops", True):
        _trim_terminal_missing(fixer)
    else:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)  # strip metals + crystal waters -> apo
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # Break the crystal disulfide BEFORE protonating, else addMissingHydrogens leaves
    # Cys57/Cys146 bonded and unprotonated -> the OXIDIZED form (see module docstring).
    strip_disulfide_bonds(fixer.topology)
    fixer.addMissingHydrogens(scfg["ph"])     # CYS -> SH (free thiols) at this pH
    assert_disulfide_reduced(fixer.topology)  # verify; never assume

    if solvate:
        forcefield = ForceField("amber14-all.xml", f"amber14/{scfg['water_model']}.xml")
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(
            forcefield,
            model=scfg["water_model"],
            padding=scfg["solvent_padding_nm"] * unit.nanometer,
            ionicStrength=scfg["ion_conc_M"] * unit.molar,
            neutralize=True,
        )
        top, pos = modeller.topology, modeller.positions
    else:
        # WT/unsolvated: Perses solvates the hybrid itself (FEP folded leg).
        top, pos = fixer.topology, fixer.positions

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        PDBFile.writeFile(top, pos, f)
    return out_path


def prepare_wt_apo(cfg: dict, variant: str, out_path: str | Path) -> Path:
    """WT (unmutated), apo/reduced, UNSOLVATED protein for the Perses FEP folded leg.

    Same chain selection + loop rebuild + apo/reduced protonation as the mutant prep, but
    the point mutation and solvation are left to Perses (see :mod:`src.fep.window`). The
    variant is still required so the correct oligomeric state / chains are selected.
    """
    return prepare_variant(cfg, variant, out_path, apply_mutation=False, solvate=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare an apo SOD1 system (Stage 1).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--wt", action="store_true",
                        help="WT, unsolvated apo protein for the Perses FEP folded leg "
                             "(no mutation, no solvent).")
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.wt:
        out = prepare_wt_apo(cfg, args.variant, args.out)
        print(f"Prepared WT apo (unsolvated) for {args.variant} -> {out}")
    else:
        out = prepare_variant(cfg, args.variant, args.out)
        print(f"Prepared {args.variant} -> {out}")


if __name__ == "__main__":
    main()
