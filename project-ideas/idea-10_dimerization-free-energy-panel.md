# Idea 10 — Dimerization Free-Energy Panel

## Research question

**Do ALS-causing mutations destabilize the SOD1 dimer interface, and can alchemical FEP
quantify that?** The parent repo computes ΔΔG of *folding* for the apo monomer. But SOD1 is a
homodimer, and the Lindberg 2005 "class 2" mutations act specifically at the **dimer
interface** (weakening the interface without destabilizing the monomer). This project computes
**ΔΔG of dimerization** — the free-energy change of assembling two monomers — for a panel of
interface and non-interface variants, and tests whether interface-localized mutations are
preferentially detected.

This is the natural **"next physics"** extension of the parent project: same engine, same
discipline, one extra leg (the dimer→two-monomer dissociation), and a distinct experimental
anchor (Hsueh et al. 2022 and Wells et al. 2021 both report SOD1 dimer ΔΔG that agree with
ITC/SEC experiments).

## Why this needs an SGE cluster

```
tasks = variants × legs × windows × replicates × (monomer | dimer) systems
      ≈ 10 × 2 × 18 × 3 × 2 = 2,160 GPU windows
```

- The dimer leg adds an **entropically hard** term (translational/rotational entropy of
  dissociation) that needs careful sampling and **many windows** — exactly the compute the
  cluster provides.
- Every window is an independent GPU job (`qsub -t`), identical in shape to the parent repo's
  Stage 3. The marginal cost over the parent project is a factor of ~2–4 in window count.
- The panel × leg × window × replicate geometry *is* the parent repo's array decode
  (`array_unit: [variant, leg, window, replicate]`) — reuse `scripts/submit_array.sh`.

## Data & tools

- **Engine:** GROMACS + pmx (parent stack). Mutation in the *dimer* (both chains) and in the
  *monomer*; the dimerization ΔΔG comes from combining the folded-dimer and folded-monomer
  alchemical legs with a dissociation reference.
- **Panel:** the parent repo's gate subset + variants with known dimer effects (Lindberg
  class-2: L144F, I104F; interface vs buried-core split is the key comparison).
- **Experimental anchors:** ITC dimerization ΔΔG (Wells 2021 cites these; Hsueh 2022 for
  computed values). The parent repo's `exp_ddg_dimer` column exists but is a carried
  cross-check — here it becomes the headline target.

## Skill prerequisites

- **Advanced.** You must be comfortable with the parent repo's full FEP machinery and with
  thermodynamic cycles for association (the monomer→dimer dissociation leg is not just "run
  pmx on a bigger system" — the reference state and symmetry factor matter).
- Understand entropic contributions to binding and why they are hard.

## Cluster budget

| Parameter | Value |
|---|---|
| Variants | 10 (5 interface, 5 buried) |
| Systems per variant | dimer + monomer |
| Windows | ~18 × 2 legs × 3 reps × 2 systems |
| Per-window | ~3–6 h GPU |
| Wall-clock on 8 GPUs | **~2–4 weeks** (the most expensive idea here) |

## Milestones

1. Reproduce a **single** dimerization ΔΔG on one well-studied interface variant (validate
   against Wells' ITC numbers before scaling).
2. Nail the dissociation leg: the reference for two free monomers must be handled correctly
   (box-size / symmetry corrections). Get one variant's error bars small enough to trust.
3. Array over the panel.
4. Aggregate per-variant ΔΔG_dimerization ± error (MBAR + cycle closure, parent discipline).
5. **The science:** does ΔΔG_dimerization separate interface variants from buried-core
   variants? Correlate with experiment. Compare dimer vs monomer ΔΔG — which mutations act
   through the interface?
6. Cross-check with the parent repo's folding ΔΔG: build a combined "monomer stability +
   dimer stability" 2-D view of the panel.

## Deliverables

- ΔΔG_dimerization table (with error bars) for the panel, vs experiment.
- **Interface-vs-buried scatter** — the headline figure.
- A 2-D map (monomer ΔΔG × dimer ΔΔG) classifying variants into Lindberg classes 1 / 2 / 1+2.
- Written verdict: can FEP on a cluster reproduce dimer-interface destabilization?

## Pitfalls

- **The dissociation leg is the trap.** Two monomers in one box carry an entropy penalty
  (ideal-gas volume) that does not cancel in the cycle unless the reference is defined
  consistently. Read how Wells/Hsueh define it before running anything.
- **Symmetry.** Homodimers have a symmetry number (2); a naive ΔΔG has a −kT·ln2 error.
  Handle it explicitly and document it.
- **Charge neutrality applies to dimer legs too** (parent config:161-166) — keep the panel
  charge-neutral, or accept and report the finite-size caveat.
- **Force-field sensitivity is doubled** (two chains, interface contacts). This is the best
  candidate for a consensus-FF treatment (Idea 02) rather than a single FF.
- **Don't launch the full array on a broken leg.** One validated interface variant end-to-end
  is worth more than 2,000 windows of a protocol with a sign error.
