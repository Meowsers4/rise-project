# Idea 06 — Binding Free-Energy Consensus (MM/PBSA + TI/MBAR)

## Research question

**Which free-energy method best ranks a small set of ligands against a protein target?**
Compute relative binding free energies for a handful of well-studied protein–ligand complexes
(e.g., inhibitors of a kinase, or known binders of a disease target) with two independent
method families:

1. **Alchemical TI/MBAR** (the rigorous, expensive way — the parent repo's GROMACS + pmx stack,
   or the `gmx bar` route), and
2. **MM/PBSA** (the cheap end-state way, `gmx_MMPBSA`).

Then compare both against **experimental binding data** and against each other. The science is
a controlled methods benchmark: does the cheap method agree with the rigorous one? Which one
is closer to experiment for this target class?

This is a natural **"upgrade path"** from the parent repo: the same `qsub -t` fan-out, the
same GROMACS machinery, but the perturbation is a ligand instead of a protein residue.

## Why this needs an SGE cluster

```
tasks = ligands × replicates × force_fields  (for TI/MBAR: + legs × windows)
```

- TI/MBAR per (ligand, FF): ~2 legs × 18 windows × 3 ns ≈ 54 windows of GPU mdrun ≈ **weeks
  single-GPU**. On 8 GPUs with a few ligands it becomes a weekend.
- MM/PBSA is cheaper (~10–20 CPU-h per complex) but still × ligands × reps → an array.
- The two method families have **different error structures**; the cluster lets you run both
  on the same ligands so the comparison is apples-to-apples.
- You can afford **multiple replicates and force fields** — which is exactly what converts a
  "docking score" into a "free energy with an error bar".

## Data & tools

- **Benchmark sets (free):** the **PDBbind refined set** (registered access, gold standard),
  or a small curated set of 5–10 protein–ligand complexes with measured Ki/IC50 (from ChEMBL
  or a paper's table).
- **Structures:** complex PDBs (bound ligand → place the ligand), ligand parameters via
  `antechamber` (AMBER) or `acpype`/`CGenFF` for GROMACS.
- **Tools:**
  - Alchemical: GROMACS `gmx bar`/`gmx ti` or **pmx** (this repo) adapted for a ligand morph.
  - MM/PBSA: **`gmx_MMPBSA`** (free, well-documented, needs AmberTools).
  - Force fields: GAFF2 (ligand) + amber99sb*-ILDNP (protein), and/or a second FF for the
    consensus idea.

## Skill prerequisites

- Comfortable with GROMACS topology building (or willing to learn — `acpype`/`antechamber`
  is the hard part).
- Understands "binding free energy", "IC50", "end-state vs alchemical".
- This is the **advanced FEP track**; do Ideas 02/05 first.

## Cluster budget

| Parameter | Value |
|---|---|
| Ligands | 5–8 |
| MM/PBSA per ligand | ~10–20 CPU-h (× 2 FFs) |
| TI/MBAR per ligand (2 FFs × 3 reps × 18 windows × 2 legs) | ~200 GPU-window-hours |
| Wall-clock on 8 GPUs | **~1–2 weeks** |

## Milestones

1. Select 5–8 complexes with clean experimental binding data; fetch PDBs + ligands.
2. Build topologies for one complex (protein: `pdb2gmx`; ligand: `antechamber`→`acpype`);
   verify the ligand does not clash in the pocket (short minimization).
3. Run **MM/PBSA** on one complex end-to-end; get the number; check sign/scale vs experiment.
4. Run **one TI/MBAR** alchemical leg on one complex to shake out the protocol (reuse the
   parent repo's window/mock machinery — it already has MBAR + provenance).
5. Array: MM/PBSA over (ligand × FF × rep); TI/MBAR over (ligand × FF × rep × legs × windows).
6. Compare: per-ligand ΔG from both methods vs experiment; rank correlation across the set.
7. The money figure: predicted vs experimental ΔG (kcal/mol) with method A, B, and the
   consensus; error bars from replicates.

## Deliverables

- Binding ΔG table (both methods, with uncertainty) vs experimental ΔG/Ki.
- A rank-order comparison: do the methods agree on which ligand binds best?
- A verdict on whether MM/PBSA is a usable cheap surrogate for this target, and where it
  breaks.

## Pitfalls

- **Ligand parameterization is the whole game.** An off-center partial charge on one atom can
  shift ΔG by 1 kcal/mol. Validate the ligand parameters (gas-phase minimization, reproduce a
  known bond length) before anything else.
- **Charge vs alchemical consistency.** If you mix AMBER (GAFF2 ligand) and GROMACS force
  fields, the charge model must match; use `acpype` consistently.
- **MM/PBSA is entropically weak** — it often gets rank order right and absolute ΔG wrong.
  Present ranks as the headline, absolute values as "estimates".
- **Don't average methods blindly.** MM/PBSA and TI are not the same estimator; averaging
  hides that. Show both plus a *documented* consensus rule.
- **Provenance (parent-repo lesson).** Stamp every result with method + FF + replicate;
  never mix unlabelled numbers.
