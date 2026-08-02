# Idea 04 — Virtual Drug-Repurposing Screen

## Research question

Can you find **approved drugs that dock well against a disease target** — i.e., a cheap
computational screen for drug repurposing? Take a target protein (e.g., an ALS-linked protein,
a viral protein like SARS-CoV-2 main protease, or a cancer kinase) and dock a library of
FDA-approved compounds against it. Rank by docking score, cluster by chemotype, and nominate
a short list of repurposing candidates.

The honest framing matters: docking scores are **ranking tools, not free energies** — the
deliverable is a prioritized shortlist with a confidence rationale, not a "drug works" claim.

## Why this needs an SGE cluster

```
tasks = number of ligands (thousands)
```

- **1 array task = 1 ligand (or a chunk of ~20 ligands).** Each docking is 5–30 min CPU
  (AutoDock Vina) or ~1–2 min GPU (GNINA). With ~3,000 FDA drugs: **~500–1,500 CPU hours** or
  a few GPU-hours on 8 GPUs.
- Docking is the canonical embarrassingly-parallel workload — no communication between
  ligands, trivially resumable, and the scoring step aggregates thousands of small results.
- Bonus: you can afford **multiple docking poses and exhaustiveness** when the cluster gives
  you 30× the cores.

## Data & tools

- **Target structure:** a PDB with a defined binding site (co-crystallized ligand → pocket),
  or a ColabFold prediction (ties to Idea 03). Choose something with a known active-site pocket.
- **Ligand library (free):**
  - **DrugBank** approved set (~2,700 compounds, free registration) — the classic repurposing
    library.
  - **ZINC** (free) subsets if you want larger/larger-varied libraries.
  - **RCSB PDB ligands** already in the pocket as a positive control.
- **Dockers:**
  - **AutoDock Vina** (free, CPU, simple `--receptor --ligand --exhaustiveness`).
  - **GNINA** (free, GPU-accelerated Vina fork — the cluster GPUs are the point).
- **Prep tools:** `rdkit` (SMILES→3D, protonation, tautomers), `obabel` or `MGLTools`
  `prepare_ligand4.py`.

## Skill prerequisites

- Python + `rdkit` to parse SMILES and slice the library.
- Basic chemistry vocabulary (SMILES, protonation, torsion).
- **No MD.** Intermediate beginner project.

## Cluster budget

| Parameter | Value |
|---|---|
| Ligands | ~3,000 (DrugBank approved) |
| Per-ligand | ~10 min CPU (Vina) or ~1–2 min GPU (GNINA) |
| Array size | 1 per ligand (or per 20-ligand chunk) |
| Wall-clock on 8 GPUs (GNINA) | **~6–12 h** |

## Milestones

1. Pick the target + pocket; download structure; extract pocket residues from the
   co-crystallized ligand.
2. Build the ligand library as SMILES → 3D SDF with `rdkit`. Add **positive controls**
   (known binders, e.g., the co-crystal ligand itself and known inhibitors) to calibrate the
   score cutoff.
3. Prepare the receptor (remove waters/heteroatoms, add H). Test **one** ligand
   interactively; sanity-check that your positive control scores well.
4. Array: `SGE_TASK_ID` → ligand file → dock → write per-ligand JSON (best score + poses).
5. Aggregate: rank ligands by score; plot score histogram; mark the positive-control score as
   the "actionable threshold".
6. **Cluster hits by 2D scaffold** (rdkit `Bemis-Murcko`) so the shortlist is chemotype
   diversity, not 50 copies of the same scaffold.
7. (Stretch) re-dock the top 50 with higher exhaustiveness + GNINA's affinity model; cross-check
   with a second docker (e.g., Smina or DOCK) as a consensus.

## Deliverables

- **Ranked repurposing table** (drug name, target score, scaffold class, known indication).
- Score-distribution histogram with the positive-control cutoff marked.
- 3D pocket–ligand pose figures for the top candidates.
- A one-paragraph "why these drugs, why this target" rationale per top hit.

## Pitfalls

- **Docking scores are not binding free energies.** Say this explicitly in every figure and in
  the write-up; a judge who does not see you know this will assume you don't.
- **Protonation/tautomer sloppiness** silently changes results. Use consistent `rdkit` rules;
  document them.
- **The pocket defines the answer.** Picking a big flat pocket will rank anything. Choose a
  well-defined site and say which one you chose.
- **Charged ligands / flexible targets.** Consider adding flexible-residue docking only for the
  top hits; full flexibility for 3,000 ligands is not worth the cost.
- **Reproducibility:** fix the random seed / exhaustiveness and record it, or two runs give two
  different top-10 lists.
