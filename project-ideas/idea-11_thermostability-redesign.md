# Idea 11 — Thermostability Redesign

## Research question

**Can consensus alchemical FEP design a set of mutations that *stabilize* a protein?**
Directed-evolution and protein-engineering projects constantly need thermostable variants of
useful proteins (industrial enzymes, biotherapeutics). This project flips the parent repo's
question from "which disease variants destabilize?" to "which designed mutations stabilize?"

Concretely: take an enzyme with a known but modest thermostability, generate a library of
candidate single- and double-point mutations (evolutionarily-informed: positions that vary in
homologs are good candidates), compute ΔΔG with the consensus-FF FEP stack, and return a
ranked list of stabilizing mutations. Validate the ranking on any mutations with published
experimental ΔTm.

## Why this needs an SGE cluster

```
tasks = candidates × legs × windows × replicates × FFs
      ≈ 20 × 2 × 18 × 3 × 2 = 4,320 GPU windows
```

- A single mutation's ΔΔG by FEP is ~2–6 GPU-days of windows. Doing 20 candidates needs the
  cluster; doing the *consensus* (2 FFs, 3 reps) makes it 2–3× more — but that is exactly the
  accuracy needed to trust a +1 kcal/mol stabilization over noise.
- The fan-out is the parent repo's Stage 3 geometry plus an FF and candidate dimension.
- The prescreen (FoldX/Idea 01) is CPU-parallel too: use it to *pre-rank* candidates cheaply
  and only FEP the top ones — a two-stage design the cluster makes affordable.

## Data & tools

- **Target protein:** an enzyme with: a PDB structure, at least one known
  thermostabilizing mutation in the literature (validation!), and a clear engineering interest
  story. Good beginner options: T4 lysozyme, barnase, a thermophilic vs mesophilic pair
  (e.g., comparing to see which residues differ).
- **Candidate generation:** B-factors / flexibility, conservation (ConSurf), FoldX scan
  (Idea 01) ranks, or published stabilizing mutations of homologs.
- **Engine:** GROMACS + pmx consensus-FF (parent stack + Idea 02's FF dimension).
- **Validation:** published ΔTm / ΔΔG for the target (e.g., T4 lysozyme has dozens of
  characterized mutations).

## Skill prerequisites

- **Advanced.** Full FEP pipeline competence (the parent repo is the prerequisite), plus
  enough protein-engineering vocabulary to talk to a mentor about thermostabilization.

## Cluster budget

| Parameter | Value |
|---|---|
| Candidates | 20 (15 single + 5 double) |
| Per-candidate | 2 FFs × 2 legs × 18 windows × 3 reps |
| Per-window | ~3–6 h GPU |
| Wall-clock on 8 GPUs | **~2 weeks** (or 1 week for singles only) |

## Milestones

1. Lock the target + a list of 20 candidates with a documented rationale each.
2. Pre-screen all candidates with FoldX (cheap, CPU) — sanity-check sign convention; flag any
   that FoldX already calls clearly destabilizing.
3. Validate the protocol on **one known-stabilizing literature mutation** — it must come out
   stabilizing in FEP, or the protocol is wrong.
4. Run the FEP array (candidates × FF × legs × windows × reps); MBAR per (candidate, FF);
   consensus per candidate.
5. Rank candidates; compare with the FoldX pre-screen; produce a top-5 shortlist with error
   bars.
6. **The science:** does the computed ranking match the literature mutation's reported
   stabilization? How many candidates pass a +1 kcal/mol threshold with error bars excluding
   zero?

## Deliverables

- **Ranked stabilizing-mutation shortlist** with ΔΔG ± error (consensus FF).
- FEP-vs-FoldX candidate agreement plot.
- Validation point: the known literature mutation vs your prediction.
- A structural rationale for the top hits (what the mutation changes at the atomic level).

## Pitfalls

- **Consensus or nothing.** Single-FF FEP has ~±1 kcal/mol noise; a "stabilizing" call from
  one FF is frequently a coin flip. This project *requires* the 2-FF consensus to be credible.
- **Double mutations are not additive.** A double mutant may stabilize more or less than the
  sum of its parts (coupling). Treat additivity as a hypothesis to *test*, not an assumption.
- **The experimental anchor is precious.** If the literature mutation doesn't come out
  stabilizing in your protocol, stop and fix the protocol before scoring the rest — otherwise
  the entire ranking inherits the error.
- **ΔTm ≠ ΔΔG.** Published thermostability is often given as ΔTm; converting Tm shifts to
  ΔΔG is model-dependent. Compare in rank space primarily.
- **Reproducibility:** keep the exact candidate list + protocol recorded before the array
  (pre-register), or the "top 5" is just whatever survived your fiddling.
