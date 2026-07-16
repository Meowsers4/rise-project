# SOD1 Variant Stability Pipeline

A GPU-cluster pipeline that computes protein-stability changes (ΔΔG of folding)
for a panel of clinically observed **SOD1** variants using rigorous **alchemical
free-energy (FEP/TI)** calculations, validated against experimentally measured
controls, then extended to interpret uncharacterized ALS variants.

> **Status:** specification / scaffolding. No code written yet. This README is the
> design contract. See [Open decisions](#open-decisions) — resolve those before
> generating implementation code.

---

## 1. What this project is (and is not)

**Goal:** produce a validated ΔΔG map across the SOD1 variant landscape that flags
which uncharacterized variants are likely destabilizing (and therefore plausibly
pathogenic) and which are not.

This is a **hybrid methods/target** project:
- *Target:* SOD1 (Cu/Zn superoxide dismutase, UniProt **P00441**), a 153-residue
  soluble homodimer implicated in familial ALS.
- *Method question:* can FEP-on-a-cluster reliably triage clinical variants when
  benchmarked against known experimental stabilities?

**In scope:** structure prep, empirical prescreen, alchemical ΔΔG, mechanism MD,
validation, variant classification, cluster orchestration.

**Out of scope (for v1):** wet-lab work, holo/metal-bound simulations (see
[§4.1](#41-apo-first-decision)), ligand docking, ML surrogate models. These are
possible later arms, not the first deliverable.

---

## 2. Core design decisions (do not silently override)

These two decisions propagate through the entire pipeline. An agent modifying this
project must treat changing either one as a scope change, not a refactor.

### 2.1 Apo-first
Simulate the **apo, disulfide-reduced** form in v1.
- **Why:** avoids metalloprotein force-field parameterization (Cu/Zn coordination),
  AND lands on the aggregation-prone, disease-relevant species. Rare case where the
  simplification is also the more biologically correct choice.
- Holo/metal-bound is a **separate parameterized arm**, not a config flag.

### 2.2 Control validation is a GATE, not a final step
Before spending cluster time on uncharacterized variants, reproduce the **known
experimental ΔΔG** of the characterized controls (A4V, G93A, G37R, …).
- If the control correlation is poor, extending to novel variants produces
  confident nonsense. Treat this as **go/no-go**.
- The control-reproduction result is itself a presentable outcome.

---

## 3. Pipeline overview

```
Stage 0  Build variant panel        (CPU, data)     controls + uncharacterized
Stage 1  Prepare structure          (CPU)           apo, protonate, solvate
Stage 2  Cheap prescreen            (CPU)           FoldX / Rosetta, ALL variants
Stage 3  Alchemical FEP  ★          (GPU, cluster)  ΔΔG — the workhorse
Stage 4  Unbiased MD                (GPU)           mechanism, subset only
Stage 5  Validate + interpret       (CPU)           controls gate, then novel
```

★ Stage 3 is the parallel cluster stage. Its job count fans out as:

```
jobs  =  variants × 2 legs × λ-windows × replicates
      ≈  40 × 2 × 18 × 5  ≈  7,200 independent GPU jobs
```

Every job is embarrassingly parallel (no inter-job communication), so throughput
scales linearly with GPU count.

---

## 4. Stages in detail

### Stage 0 — Build and stratify the variant panel
- Sources: **ClinVar**, **ALSoD** (ALS-specific), cross-ref **gnomAD** for frequency.
- **Numbering caution:** mature SOD1 has the initiator Met removed → 153 residues.
  Literature variant names (A4V, G93A) use mature numbering. Get the offset right
  once, centrally, or every downstream mapping breaks.
- Stratify into four buckets:
  1. **Positive controls** — experimentally measured ΔΔG / ΔTm (non-negotiable).
  2. **Negative controls** — common/benign, should not destabilize.
  3. **Pathogenic-but-uncharacterized** — clinical label, no structural work.
  4. **VUS** — variants of uncertain significance (the payoff targets).
- Output: `data/variants.csv` with columns
  `variant,mature_pos,wt_aa,mut_aa,bucket,exp_ddg,exp_source,clinvar_id`.

### Stage 1 — Prepare structure
- Start from a high-resolution human SOD1 crystal structure (PDB).
- **Apo** (see §2.1). Strip metals; model disulfide-reduced Cys57/Cys146.
- Monomer vs dimer: **dimer** for interface-adjacent variants, **monomer** for
  buried-core variants — decided per variant, recorded in `variants.csv`.
- Standard prep: fix missing atoms/loops, assign protonation (PROPKA or H++),
  solvate, neutralize, add ions.
- Deliver as a **parameterized script**: `variant_id -> prepared_system`. It runs
  dozens of times; no manual steps.

### Stage 2 — Cheap empirical prescreen (CPU, runs on everything)
- **FoldX BuildModel** and/or **Rosetta cartesian_ddg** on ALL variants.
- Three jobs: first ΔΔG estimate, prioritization ranking, and an independent
  cross-check against FEP later.
- Validate these against controls too — if FoldX can't reproduce known cases,
  you learn it cheaply.

### Stage 3 — Alchemical ΔΔG (the cluster workhorse) ★
- Thermodynamic cycle: run the alchemical mutation (X→Y) in the **folded protein**
  and in an **unfolded reference** (solvated capped tripeptide). Difference of the
  two legs = ΔΔG_folding.
- Each leg → **12–24 λ-windows**, softcore potentials for appearing/disappearing
  atoms, **≥3 (ideally 5) replicates** per window with distinct seeds.
- Replicate discipline is where people cut corners; error bars are mandatory here.
- **Framework — one of (see [Open decisions](#open-decisions)):**
  - OpenMM + Perses
  - GROMACS + pmx
  - AMBER TI
- Free-energy estimator: **MBAR**. Check convergence (forward/backward hysteresis,
  cycle closure).

### Stage 4 — Unbiased MD for mechanism (complementary)
- FEP gives a number; MD gives the story. Run longer replicate plain-MD on the
  subset FEP flags as interesting.
- Analyze: per-residue **RMSF**, contact-map / salt-bridge loss, **SASA** change,
  dimer-interface integrity, local-unfolding proxies.
- **Caveat:** do not expect full unfolding in accessible timescales — read proxies
  and equilibrium shifts, not a complete denaturation event.

### Stage 5 — Validation and interpretation
- **GATE:** correlate computed FEP ΔΔG vs experimental controls (Pearson/Spearman,
  low cycle-closure hysteresis). Go/no-go.
- Cross-check FEP vs FoldX/Rosetta for consistency.
- Then classify uncharacterized variants; ask whether computed destabilization
  separates clinically pathogenic from benign labels.
- Output: `results/ddg_map.csv` + validation plots.

---

## 5. Cluster orchestration

- **Workflow manager:** Snakemake or Nextflow — the DAG
  (prep → prescreen → FEP windows → analysis) must be reproducible and restartable.
- **Scheduler:** SGE/OGS (BU SCC) **job arrays** (`qsub -t`) — one array task per
  `(variant, leg, window, replicate)`. GPUs requested by compute capability
  (`-l gpus=1 -l gpu_c=7.0`), never by model; the scheduler assigns the device via
  `CUDA_VISIBLE_DEVICES`, so the FEP code must not hardcode it.
- **Checkpoint everything** — windows die; never restart from zero.
- **Trajectory retention policy up front** — full FEP trajectories are large.
  Keep free-energy estimates + analysis outputs; downsample or discard raw frames.

---

## 6. Proposed repository layout

```
sod1-fep/
├── README.md                  # this file
├── env/
│   ├── environment.yml        # conda: mdtraj, openmm/gromacs, pymbar, etc.
│   └── modules.md             # cluster module-load notes
├── data/
│   ├── variants.csv           # Stage 0 output (the source of truth)
│   └── structures/            # prepared PDBs, per variant
├── config/
│   └── pipeline.yaml          # global params incl. cluster block (SGE): project, gpu_c, walltime
├── workflow/
│   └── Snakefile              # or main.nf
├── src/
│   ├── panel/                 # Stage 0: pull + stratify variants
│   ├── prep/                  # Stage 1: variant_id -> prepared_system
│   ├── prescreen/             # Stage 2: FoldX / Rosetta wrappers
│   ├── fep/                   # Stage 3: setup, run, MBAR analysis
│   ├── md/                    # Stage 4: unbiased MD + analysis
│   └── analysis/              # Stage 5: validation, plots, classification
├── scripts/
│   └── submit_array.sh        # SGE array submission helper (qsub -t)
├── results/
│   ├── ddg_map.csv
│   └── figures/
└── tests/
    └── test_numbering.py      # guard the mature-numbering offset
```

---

## 7. Suggested tech stack

| Concern            | Choice                                              |
|--------------------|-----------------------------------------------------|
| Language           | Python 3.11+                                         |
| MD / FEP engine    | OpenMM+Perses **or** GROMACS+pmx **or** AMBER TI     |
| Empirical ΔΔG      | FoldX, Rosetta cartesian_ddg                         |
| Free-energy math   | pymbar (MBAR)                                        |
| Trajectory analysis| MDTraj / MDAnalysis                                  |
| Workflow           | Snakemake or Nextflow                                |
| Scheduler          | SGE/OGS (BU SCC), job arrays (`qsub -t`)              |
| Structure prep     | PDBFixer, PROPKA/H++                                 |

---

## 8. Milestones

1. **M0** — variant panel built and stratified; numbering test passes.
2. **M1** — structure-prep script runs end-to-end on one control variant.
3. **M2** — empirical prescreen over all variants; ranked table produced.
4. **M3** — FEP runs for the control set; **validation gate** evaluated.
   *(No further compute until this gate passes.)*
5. **M4** — FEP extended to uncharacterized variants; ΔΔG map produced.
6. **M5** — mechanism MD on flagged subset; final report + figures.

---

## 9. Open decisions

An agent should ask the user to resolve these before writing implementation code:

- **FEP framework:** OpenMM+Perses vs GROMACS+pmx vs AMBER TI? (Shapes all of Stage 3
  and the SGE array structure.)
- **GPU count / partition names / walltime limits** on the target cluster? (Sizes the
  panel and the array chunking.)
- **Panel size for v1:** how many variants (≈30–50 suggested)?
- **Monomer vs dimer default**, and the per-variant rule for choosing.
- **Starting PDB** (which structure, which resolution).
- **λ-window count (12/18/24)** and **replicate count (3/5)** — cost vs precision.
- **Alternative target fallback:** if the SOD1 metal chemistry proves troublesome
  even in apo form, **transthyretin (TTR)** offers the same panel+FEP design with no
  metals (tetramer stability is the disease-relevant quantity). Keep as a documented
  Plan B.

---

## 10. References to gather (for the agent to populate)

- SOD1 experimental stability data for the control variants (literature ΔTm/ΔΔG).
- ClinVar + ALSoD variant exports.
- Chosen FEP framework's tutorial / protocol paper.
- pmx or Perses hybrid-topology setup docs.
