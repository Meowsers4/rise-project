# Idea 27 — Co-Evolution Coupling vs FEP: Do Evolutionary Signals Predict Physics?

## Research question

Evolutionary sequence analysis (e.g., **EVcouplings / DCA**: direct coupling analysis on a
multiple sequence alignment) identifies residue pairs that must co-vary to preserve
structure/function — a purely sequence-based, "free" signal. FEP computes the physical free
energy of mutations — an expensive, first-principles signal. The question that almost nobody
tests cleanly:

**Do co-evolutionary couplings predict which mutations destabilize a protein, and how does
that compare to FEP (or FoldX) on the same variants?**

The interesting part is the *disagreement*: if DCA says a residue is critically coupled but
FEP says the mutation is neutral, either (a) FEP's model misses an allosteric/kinetic effect,
or (b) DCA is detecting functional (not stability) constraint. Resolving that tension is the
science.

## Why this needs an SGE cluster

```
tasks = FEP windows for a variant panel (expensive) + DCA is cheap (CPU) + a scan of pairs
```

- The FEP side is the parent repo's array: a panel of variants (say 20–40), each a
  variant × leg × window × replicate fan-out. That's the expensive half and only the cluster
  can do it at a meaningful panel size.
- DCA is cheap but needs a deep MSA (hundreds of sequences) and produces a *full coupling
  matrix* over ~N² residue pairs — computing per-pair coupling statistics and matching them to
  mutation positions is parallelizable (array over pairs).
- You also run a **FoldX scan** (cheap) over the same panel so you can compare all three
  signals (DCA, FoldX, FEP) on the same variants.

## Data & tools

- **DCA/EVcouplings:** the EVcouplings server or `plmc`/`ccmgen` tools (free); needs a deep
  MSA (UniProt blast, or a curated MSA for SOD1 — plenty exist). Output: per-pair coupling
  scores, per-residue conservation, predicted mutational effects.
- **FEP:** parent repo engine, GROMACS + pmx.
- **Cheap physics:** FoldX scan (Idea 01).
- **Validation:** experimental ΔΔG (variants.csv; a DMS).

## Skill prerequisites

- Comfortable with the parent repo's FEP machinery.
- Basic understanding of multiple sequence alignments and what a coupling score means.
- Intermediate-advanced (runs the full spectrum cheap→expensive).

## Cluster budget

| Parameter | Value |
|---|---|
| Variant panel | 20–40 |
| FEP windows | ~20–40 × 108 = 2,200–4,300 |
| FoldX scan | ~3,000 (cheap, shares Idea-01) |
| DCA | CPU minutes once MSA is built |
| Wall-clock on 8 GPUs | **~2–3 weeks** (the FEP is the cost) |

## Milestones

1. Build a deep MSA for the target; run DCA → per-pair couplings and per-residue
   conservation.
2. Define the panel: mutations at (a) high-coupling, (b) low-coupling positions.
3. Run the cheap signals: FoldX scan + DCA's predicted effects.
4. Run the FEP panel (parent array); MBAR per variant.
5. **The comparison:** do DCA couplings rank variant severity like FEP does? Compute
   correlation between DCA score, FoldX ΔΔG, FEP ΔΔG, and experiment (where known).
6. **The tension analysis:** for variants where DCA and FEP disagree, look structurally —
   what distinguishes them (interface? allosteric path? functional site?).
7. Verdict: is DCA a free pre-screen for FEP, and when does it mislead?

## Deliverables

- **Three-signal comparison** (DCA vs FoldX vs FEP, against experiment) — the money figure.
- A coupling-vs-ΔΔG scatter with the disagreement cases highlighted.
- Structural case studies of 3–5 disagreement variants (the most interesting scientific
  content).
- A recommendation: can DCA cheaply rank the *uncharacterized* VUS before FEP (a pipeline
  decision for the parent repo).

## Pitfalls

- **DCA measures functional/evolutionary constraint, not just stability.** Strong coupling
  can be about function (e.g., an active site) rather than folding; that's WHY disagreement is
  interesting, but don't call a DCA hit a "stability prediction."
- **MSA depth matters.** Shallow MSAs give garbage couplings; check the effective sequence
  count and report it.
- **Not every coupled pair is a stability coupling.** Many DCA top pairs are in contact
  (structural), but some reflect function/allostery — annotate pair types (contact vs distant).
- **Charge and phylogeny confounds.** Phylogenetic non-independence can fake couplings;
  consider reweighting the MSA (DCA tools do this; document it).
- **FEP cost is real.** A 40-variant panel is 3+ weeks on 8 GPUs; start with 20 and prioritize
  the disagreement-testing variants (Idea 17's active learning is a natural pairing).
