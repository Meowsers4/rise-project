# Idea 01 — In-Silico Saturation Mutagenesis Scan

## Research question

Map the **complete stability landscape** of a single protein — every single-point mutation at
every position (~19 × N mutations) — with an empirical free-energy function. Which positions
are mutation-intolerant (stability "hotspots")? Does the in-silico landscape match
evolutionary conservation, and can we identify positions that are unusually tolerant (a
"mutability hotspot") where most substitutions are neutral?

A strong framing: **"Predicting protein evolvability from first principles"** — proteins that
tolerate many mutations can explore more sequence space, which is how evolution and directed
evolution work. Compare your computed landscape against conservation scores (from a multiple
sequence alignment) and, if available, a published experimental deep mutational scan.

## Why this needs an SGE cluster

The fan-out is trivial and perfect for arrays:

```
tasks = N_residues × 19 = 153 × 19 ≈ 2,900  (SOD1), or ~4,000 for a 200-res protein
```

- **1 array task = 1 mutation.** No shared state, no inter-task communication, no memory
  pressure. This is the textbook use case for `qsub -t 1-3000`.
- FoldX `BuildModel` takes ~2–5 min per mutation single-threaded → **~2,900 × 4 min ≈ 200 CPU
  hours**. On a laptop that is a week of continuous compute; on the cluster it is a few hours
  of wall-clock with `-pe omp 4` and 30+ concurrent tasks.
- The result is a **2-D matrix (position × amino acid)** that only exists because you computed
  thousands of independent points. No single-machine shortcut gives you the full landscape.

## Data & tools

- **Protein:** any well-resolved monomeric structure. Good beginner choices: SOD1 (3ECU —
  already in this repo!), T4 lysozyme, barnase (the FEP benchmark standard), GFP, or a
  therapeutic antibody domain.
- **Tools:**
  - **FoldX** (`BuildModel`) — fast empirical ΔΔG. Academic license required; the BU SCC
    usually has a licensed install — record the path in config, don't hardcode.
  - **Rosetta `cartesian_ddg`** — slower, more accurate; good as a cross-check on a subset.
  - **MutateX** (Tiberti et al. 2022, PMID 35323860) automates the FoldX scan and ensemble
    averaging — strongly recommended to avoid writing this plumbing yourself.
- **Conservation:** generate an MSA with `blastp` against nr, or use precomputed ConSurf /
  EVcoupling / UniProt conservation data.
- **Validation:** published experimental DMS datasets (e.g., Fowler & Fields 2014 dataset
  collection, or specific proteins like GBP1 or TEM-1 β-lactamase).

## Skill prerequisites

- Python (pandas, matplotlib) to orchestrate and plot.
- Basic structural biology vocabulary (residue, side chain, buried vs surface).
- **No** MD/FEP knowledge required — this is the gentlest on-ramp.

## Cluster budget

| Parameter | Value |
|---|---|
| Mutations | ~2,900–4,000 |
| Per-task time | ~2–5 min CPU |
| Array size | 1 per mutation |
| Wall-clock (30 concurrent, `-pe omp 4`) | **~6–12 hours** |

## Milestones

1. Pick the protein; download the PDB; verify it has no missing residues at your positions
   (or let FoldX repair).
2. Get FoldX running interactively on ONE mutation; write the output JSON you will collect.
3. Build the array script: decode `SGE_TASK_ID` → (position, amino acid); run `BuildModel`
   (optionally 5 runs and average, `prescreen.runs_per_variant` style); write
   `results/mutscan/{pos}_{aa}.json` with provenance.
4. `qsub -t 1-<N>`; poll `qstat`; collect outputs.
5. Assemble the ΔΔG matrix; make a heatmap (position × amino acid).
6. Correlate with conservation: are intolerant positions conserved? Compute per-position mean
   |ΔΔG| vs MSA entropy.
7. (Stretch) compare a subset to Rosetta `cartesian_ddg` and/or experimental DMS.

## Deliverables

- **Stability heatmap** for the full protein (the money figure).
- Scatter/rank plot: in-silico ΔΔG vs conservation entropy.
- A ranked list of the most intolerant and most tolerant positions, with a structural
  explanation (buried core vs surface loop) for each.
- A short "evolvability map" interpretation: which regions can evolution explore freely?

## Pitfalls

- **Numbering.** If you reuse this repo, remember SOD1 uses **mature numbering** (offset in
  config). For a new protein, record which chain/residue-numbering scheme you used in the
  results, or a judge will get confused.
- **Repair loops.** FoldX needs a complete, repaired structure. Run its repair on the crystal
  structure before the scan, or every mutation is evaluated on a slightly-broken protein.
- **Sign convention.** Positive ΔΔG = destabilizing (matching `validation.ddg_sign`). State it
  once in the README of your results.
- **Parallel pdb-name collisions.** FoldX writes files into its working directory; give each
  task its own scratch dir (e.g., `${TMPDIR}` or a per-task folder) or concurrent tasks will
  clobber each other.
- **Don't trust a single run.** Average ≥3–5 BuildModel runs per mutation (FoldX is noisy,
  ~0.5 kcal/mol) — this is exactly what `prescreen.runs_per_variant` does in the parent repo.
