"""Stage 1 -- prepare an apo, disulfide-reduced SOD1 system for one variant.

A parameterized ``variant_id -> prepared_system`` step (README §4 Stage 1), run once
per variant with no manual steps. Pipeline (PDBFixer + OpenMM ``Modeller``):

1. Load ``structure.starting_pdb`` (3ECU) and keep only the chains for this variant's
   oligomeric state (monomer -> ref chain; dimer -> ref + geometric partner).
2. Apply the point mutation in every kept chain (a homodimer carries it in both).
3. Rebuild 3ECU's disordered metal-binding loops (D5) via PDBFixer missing-residue/
   atom addition -- internal gaps only, not terminal tails.
4. Strip all heteroatoms (metals + crystal waters) -> apo (§2.1).
5. Protonate at ``structure.ph``; cysteines get HG (free thiols) -> disulfide REDUCED.
6. Solvate (``water_model``) with ``solvent_padding_nm`` padding and neutralize to
   ``ion_conc_M``.

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


def _trim_terminal_missing(fixer) -> None:
    """Drop terminal missing residues so we rebuild internal loops only, not new tails."""
    chains = list(fixer.topology.chains())
    for key in list(fixer.missingResidues.keys()):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]


def prepare_variant(cfg: dict, variant: str, out_path: str | Path) -> Path:
    """Build the apo, reduced, solvated system for ``variant`` and write it to PDB."""
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

    # mutate every kept chain (homodimer carries the mutation in both subunits)
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
    fixer.addMissingHydrogens(scfg["ph"])     # CYS -> SH (reduced) at this pH

    forcefield = ForceField("amber14-all.xml", f"amber14/{scfg['water_model']}.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        forcefield,
        model=scfg["water_model"],
        padding=scfg["solvent_padding_nm"] * unit.nanometer,
        ionicStrength=scfg["ion_conc_M"] * unit.molar,
        neutralize=True,
    )

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare an apo SOD1 system (Stage 1).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = prepare_variant(cfg, args.variant, args.out)
    print(f"Prepared {args.variant} -> {out}")


if __name__ == "__main__":
    main()
