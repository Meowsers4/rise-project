# Idea 07 — Peptide-Library Interface Scan

## Research question

**Which short peptides bind best to a therapeutically relevant protein interface?** Take a
target interface (e.g., a protein–protein interaction implicated in disease — a homodimer
interface like SOD1's, a PDZ domain, an MHC-peptide groove, or a viral entry receptor like
ACE2) and dock a **systematic peptide library** against it. Rank by docking score, then
re-rank the top candidates with short MD or MM/PBSA.

A strong framing: **"In-silico screening of a peptide library for competitive binders at a
protein–protein interface."** Peptides are attractive as therapeutics and this is a project
whose entire scale comes from the cluster (thousands of small docking jobs, then a handful of
expensive validation runs).

## Why this needs an SGE cluster

```
tasks = library_size  (thousands of peptides)
```

- **1 array task = 1 peptide** (or a few). Docking a 5–8 residue peptide is ~10–30 min CPU
  → 2,000 peptides ≈ **700 CPU-hours**; the cluster makes this a day.
- The two-stage design needs the cluster twice: a wide cheap screen (docking) then a narrow
  expensive screen (MD/MM-PBSA) on the top ~20 — only affordable because the wide stage ran
  in parallel.
- Peptides are flexible (many torsions), so you want multiple docking poses per peptide and
  more exhaustiveness — the cluster pays for that too.

## Data & tools

- **Target:** any well-defined interface. Good candidates: SOD1 dimer interface (builds on
  the parent repo), the ACE2–spike interface, a PDZ domain, a SH2/SH3 domain, or an
  E3-ligase–substrate interface.
- **Library:** enumerate systematically (all di/tripeptides: 400/8,000), or a designed
  library (e.g., all single-point variants of a known peptide binder → "alanine scan in
  silico"), or a random subset of longer peptides.
- **Tools:**
  - Docking: **AutoDock Vina** or **HADDOCK** (HADDOCK excels at peptide-protein).
  - Prep: `rdkit` / `pep2Dock`-style builders; `MODELLER` for peptide templates if needed.
  - Validation: **GROMACS MD** (short 10–50 ns per top hit) or `gmx_MMPBSA`.

## Skill prerequisites

- Python, basic chemistry (residues, charge, peptide bond).
- Some experience with a docking tool (Idea 04 is a gentler intro).

## Cluster budget

| Parameter | Value |
|---|---|
| Library | 2,000–8,000 peptides |
| Per-peptide docking | ~10–30 min CPU |
| Top-20 validation (MD/MMPBSA) | ~5–20 GPU-h each |
| Wall-clock | **~2–4 days** (docking) + 1–2 days (validation) |

## Milestones

1. Choose the interface; download the structure; define the binding region (contact residues).
2. Build the peptide library as 3D PDB/SDF files (`rdkit` → peptide builders).
3. Dock ~50 peptides locally; sanity-check that a **known binder** (if the interface has one)
   scores well.
4. Array over the library; collect best-pose scores.
5. Rank; apply a cutoff; select top 20–50 for a **second, more exhaustive docking round**
   (higher exhaustiveness, flexible peptide) — still cheap on the cluster.
6. Validate the top 10–20 with short GROMACS MD (complex stability: RMSD, contact occupancy)
   or MM/PBSA.
7. Present: ranked table, pose figures, and a "structural rationale" for the top hits
   (which interface residues they engage).

## Deliverables

- **Ranked peptide library** (score, scaffold, pose).
- Validation table (top peptides re-scored by MD/MMPBSA).
- 3D interface figures showing the top peptide bound.
- A one-line "design rule" if one emerges (e.g., "pos 3 must be hydrophobic").

## Pitfalls

- **Flexibility.** Peptides are floppy; a single rigid pose from Vina can mislead. Always
  follow up top hits with an MD or higher-exhaustiveness round — that two-stage design is the
  scientific core.
- **Interface choice defines the answer.** Pick a well-formed pocket; a flat protein surface
  docks everything anywhere.
- **Scoring artifacts:** charged termini on peptides bias electrostatics. Neutralize caps
  (ACE/NME) and document it.
- **Enumerate, don't random.** A designed library (all variants of a known binder) is far more
  interpretable than a random one — you can extract a per-position preference matrix from the
  scan (mini saturation mutagenesis at a binding interface).
- **Don't overclaim binding.** Docking + 50 ns MD suggests hypotheses about binding, not
  affinity. Report ranks, not Ki.
