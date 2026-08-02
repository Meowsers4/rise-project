# Idea 03 — Batched ColabFold / AlphaFold Proteome Scan

## Research question

Predict high-confidence structures for **every protein in an understudied family or organism**
and answer a biological question with them: Do these proteins share a fold? Are there
uncharacterized domains with known structural neighbors? Can predicted structures reveal a
function (via fold-matching to proteins of known function)?

A crisp framing: **"Mining an uncharacterized protein family for folds and functions."**
Examples:
- All proteins of a small pathogen's proteome lacking a PDB structure.
- A large enzyme family (e.g., a kinase-like family with many uncharacterized members).
- All SOD1 homologs across species → do disease-relevant structural features (dimer interface
  loops, metal-binding geometry) hold in homologs?

## Why this needs an SGE cluster

```
tasks = number of proteins to predict (hundreds)
```

- **Each ColabFold prediction is one GPU task** (~15–45 min on a single GPU for a typical
  protein, more for long ones). AlphaFold/ColabFold is *not* multi-GPU-per-job; the
  parallelism is **one protein per GPU**, exactly what an array provides.
- Predicting 300 proteins back-to-back on one laptop GPU ≈ 150+ hours (and ColabFold's
  AlphaFold2 model alone is several GB — a laptop fight). On 8 cluster GPUs this is a
  **1–2 day** wall-clock run with `qsub -t 1-300`.
- The downstream analysis (fold matching, clustering) is also trivially parallelizable when
  you have hundreds of structures.

## Data & tools

- **ColabFold** (`colabfold_batch`) — free, MIT-licensed, runs AlphaFold2 via MMseqs2 MSA
  search. GPU-accelerated; the standard choice for batch prediction.
- **Sequences:** UniProt family queries, or a proteome FASTA (e.g., from UniProt or NCBI).
  For "uncharacterized" you want sequences with no experimental structure in PDB.
- **Fold matching:** Foldseek (fast), DALI (slower), or CATH/ECOD database search.
- **MSA source:** ColabFold's MMseqs2 server search is built in; no local BLAST DB needed.

## Skill prerequisites

- Python to slice FASTA files and build the sequence list.
- Basic bioinformatics vocabulary (FASTA, UniProt ID, protein family).
- **No MD knowledge.** This is a beginner-friendly GPU-array project.

## Cluster budget

| Parameter | Value |
|---|---|
| Proteins | 200–400 |
| Per-task | ~15–45 min GPU (protein-dependent) |
| Array size | 1 per protein |
| Wall-clock on 8 GPUs | **~1–2 days** |

## Milestones

1. Choose the family/proteome; download FASTA; filter to proteins without experimental
   structures and with sensible lengths (e.g., 60–600 residues to keep runtime sane).
2. Install ColabFold (conda); run **one** protein interactively on a GPU node to learn the
   flags and output format.
3. Write the array script: `SGE_TASK_ID` → one FASTA entry; run `colabfold_batch --num-recycle 3`
   into a per-task output dir; emit a status JSON.
4. Submit `qsub -t 1-<N>`; collect predicted PDBs.
5. Score quality: keep predictions with pLDDT > 70 (confident cores); note low-confidence
   regions as flexible/unstructured.
6. Fold-matching: Foldseek-search every predicted structure against a database (e.g., PDB100);
   tabulate "best structural neighbor + function" per protein.
7. Write up: what fraction got confident structures, what functions were matched, any novel
   folds / conserved structural cores.

## Deliverables

- Predicted-structure gallery (colored by pLDDT) for the family.
- A table: protein → best fold match → inferred function → confidence.
- A structural-clustering dendrogram of the family (are they one fold or many?).
- A short functional-annotation story for the previously-unknown members.

## Pitfalls

- **GPU memory.** Some ColabFold models exceed 16 GB for long proteins; keep lengths bounded
  or the task dies at the top of the array. Use `--max-msa`/recycling defaults and test the
  longest sequence first.
- **Disk.** Each prediction writes ~10s of MB (model params are shared, cached). Use
  `trajectory_retention`-style discipline: keep the PDB + pLDDT, drop the model cache.
- **Not all proteins are globular.** Expect a fraction of intrinsically-disordered or
  multi-domain failures (low pLDDT); report them rather than deleting them.
- **Homologs are not a validation set.** Matched folds suggest function; they don't prove it.
  Phrase claims as hypotheses, and note where experimental data would be needed.
