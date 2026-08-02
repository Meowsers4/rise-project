# Idea 13 — Host–Pathogen Protein-Protein Docking Scan

## Research question

**Which pathogen proteins interact with human proteins, and where?** Pathogens (viruses,
bacteria, parasites) subvert host biology by binding host proteins. With predicted structures
now abundant (Idea 03: ColabFold every pathogen protein), you can run a **systematic
protein-protein docking scan**: every pathogen protein vs a panel of human proteins (or vs the
human interactome subset), score the interactions, and rank the most plausible host-pathogen
contacts for follow-up.

Framings:
- *"Predicting novel host factors targeted by a pathogen's proteome."*
- *"Which SARS-CoV-2 proteins, beyond spike and NSPs, have human interaction partners?"*
- *"Mining an uncharacterized pathogen proteome for interaction candidates."*

This is a **hypothesis-generation** project: docking scores rank candidates; you present a
shortlist with structural rationale, not confirmed interactions.

## Why this needs an SGE cluster

```
tasks = pathogen_proteins × host_proteins  (e.g., 20 × 25 = 500)
```

- **1 array task = 1 (pathogen, host) pair.** Each protein-protein docking is ~1–3 h CPU
  (HADDOCK) — 500 pairs ≈ **1,000 CPU-hours**; the cluster makes this days.
- Protein-protein docking is far slower than small-molecule docking and benefits enormously
  from parallelism — and from being able to test *all* pairs rather than a pre-screened few.
- You can afford **multiple docking modes** (rigid-body + refinement) per pair and a
  consensus score across tools, which is the difference between noise and a credible ranking.

## Data & tools

- **Pathogen proteins:** ColabFold predictions (Idea 03) or PDB structures of the pathogen
  proteome.
- **Host proteins:** a curated set (e.g., known immune/viral-relevant human proteins, or a
  subset of the human interactome). UniProt sequences → structures (ColabFold) if no PDB.
- **Docking tools (free):**
  - **HADDOCK** (data-driven docking; best for protein-protein; free server or local install;
    needs an interface/flexibility definition).
  - **ClusPro** (free server) or **Z-Dock/ZDOCK** as a second opinion.
- **Scoring/validation:** solvation/electrostatics scores, cluster size (ClusPro), and
  (stretch) short MD or MM/PBSA on the top pairs (Idea 06 tooling).

## Skill prerequisites

- Python + basic structural biology.
- Comfort running a docking server/package and parsing scores.
- Intermediate; lighter than the FEP ideas.

## Cluster budget

| Parameter | Value |
|---|---|
| Pairs | ~500 (20 pathogen × 25 host) |
| Per-pair | ~1–3 h CPU (HADDOCK) |
| Array size | 1 per pair |
| Wall-clock (30–50 concurrent) | **~2–4 days** |

## Milestones

1. Assemble the pathogen + host protein lists with structures.
2. **Positive control:** take one known host-pathogen interaction (e.g., spike–ACE2, or a
   published pathogen–host pair) and confirm your docking protocol ranks it high. This is
   non-negotiable before the full scan.
3. Define per-pair interface residues (or run ab-initio mode) — record the choice.
4. Array over pairs; collect per-pair top-cluster scores.
5. Rank pairs; apply a score threshold calibrated on the positive control.
6. Cluster top hits by function (Gene Ontology of host partners) — does the pathogen
   preferentially target immunity, trafficking, translation?
7. (Stretch) Short MD / MM/PBSA on the top 5–10 pairs to add a dynamics-based filter.

## Deliverables

- **Ranked host-pathogen interaction table** (pathogen protein, host protein, score, cluster
  size, function).
- Positive-control calibration plot.
- Functional enrichment of the top host partners (which pathways are targeted?).
- 3D complex models for the top candidates.

## Pitfalls

- **Calibrate on a known interaction or nothing is meaningful.** Docking scores are
  tool-dependent; without a positive control you cannot even define "high score".
- **Flexibility.** Protein-protein interfaces are conformationally variable; rigid docking
  misses induced fit. Use refinement mode and treat results as candidates, not facts.
- **Scoring artifacts:** large proteins get big favorable scores; normalize per-interface-area
  or per-residue, or your ranking is "biggest proteins first".
- **Homology leakage.** If pathogen and host proteins are actually homologs (e.g., both
  globins), you're docking similar folds, not "interacting" ones. Filter obvious homologs or
  flag them.
- **These are predictions.** The deliverable is a ranked hypothesis list with a structural
  rationale — say that clearly, and note what experimental test would confirm a hit.
