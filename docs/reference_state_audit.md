# Reference-state audit: are the experimental controls comparable to what we simulate? — 2026-09-12

The last week-1 deliverable of [`post_pivot_review.md`](post_pivot_review.md) §5, and the
question [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) §7 explicitly could not reach:
reanalysis proved the arithmetic is right; it cannot prove the calculation and the experiment
describe the same thermodynamic state.

**Headline: they do not.** Every gate control was measured with the **Cys57–Cys146 disulfide
intact**. We simulate the **disulfide-reduced (2SH)** form, deliberately and by rule. The gate
has been comparing apo-**SS** experiment against apo-**2SH** calculation on all eight points.

A second finding runs the other way: the three values traceable to a primary table
(A4V, G93A, I113T) **re-derive exactly** from Lindberg 2005. The panel's arithmetic is sound.
The problem is the state, not the bookkeeping.

This is a **scope conflict for the user to resolve**, per `CLAUDE.md` ("surface the conflict
instead of silently choosing"). Nothing in `CLAUDE.md`, `README.md` or `config/pipeline.yaml`
has been changed on the strength of it. See §7.

---

## 1. The chain of custody

`data/variants.csv` cites two sources per gate variant: `Kumar2017_TableS1` plus a primary
paper. **Kumar 2017 is a compilation, not a measurement.**

| | |
|---|---|
| Paper | Kumar V, Rahman S, Choudhry H, Zamzami MA, Jamal MS, Islam A, Ahmad F, Hassan MI (2017). *Computing disease-linked SOD1 mutations: deciphering protein stability and patient-phenotype relations.* Sci Rep 7:4678. |
| DOI / PMC | [10.1038/s41598-017-04950-9](https://doi.org/10.1038/s41598-017-04950-9) · [PMC5498623](https://pmc.ncbi.nlm.nih.gov/articles/PMC5498623/) |
| What it is | A benchmark of 8 stability predictors (PoPMuSiC, FoldX, mCSM, …) against a **compiled** experimental set: **54 apo-monomer** and **33 holo-dimer** ΔΔG values. |
| Where its numbers came from | Quoted: *"These data have been taken from Vassall et al. Nordlund and Oliveberg, Lindberg et al. Stathopulos et al., and Bystrom et al."* |
| Experimental conditions stated | **None.** No pH, temperature, method, construct or redox state for the compiled values. |
| Stated uncertainty | *"the methodological error in the experimental ΔΔG is on the order of up to ~0.3 kcal/mol."* |

Two consequences.

**The uniform `exp_ddg_err: 0.3` on every panel row is Kumar's blanket figure**, not a
per-measurement uncertainty. It is one number covering five labs, two denaturants, two
temperatures and both calorimetric and chemical methods. It should not be read as the error on
any individual control, and a gate RMSE of 1.5 is not "3–5σ" in any meaningful sense.

**The compilation is where the conditions were lost.** Kumar's set is internally
heterogeneous; our panel inherited the flattened column, and `measured_state: apo_monomer` on
all eight rows is the panel's own label, not something Kumar asserts.

## 2. Per-source audit

| gate variants | primary source | construct | oligomer | metal | **disulfide** | conditions | method |
|---|---|---|---|---|---|---|---|
| A4V, G93A, I113T | Lindberg 2005 PNAS 102:9754 | C6A/C111A **+ F50E/G51E** | **monomer** | apo (10 mM EDTA) | **intact / oxidized** | 10 mM MES **pH 6.3**, **25 °C** | urea equilibrium + kinetics |
| F64A, I149A, I18V | Nordlund & Oliveberg 2006 PNAS 103:10218 | C6A/C111A **+ F50E/G51E** | **monomer** | apo (10 mM EDTA) | **intact / oxidized** | 10 mM MES **pH 6.3** | urea |
| G93S, G93V | Stathopulos 2006 JBC 281:6184 | pseudo-WT, "no free cysteines" (C6A/C111A); **no interface mutations reported** | **dimer** (see §5) | apo and holo | intact | **DSC — thermal** | differential scanning calorimetry |

### Verified quotations

Lindberg 2005, on the redox state — the decisive sentence:

> "All experiments were done under oxidizing conditions, with the intramolecular disulphide
> linkage between C57 and C146 kept intact."

and on construct and buffer:

> "the double mutant C6A/C111A (SODpwt) as pseudo-WT for further mutant analysis" …
> "the monomeric variant F50E/G51E/C6A/C111A" … "the standard buffer was 10 mM Mes (pH 6.3)
> with 10 mM EDTA to maintain the proteins metal free"

Nordlund & Oliveberg 2006 is the same laboratory, the same pWT monomer construct
(C6A/C111A/F50E/G51E) and the same 10 mM MES pH 6.3 / 10 mM EDTA buffer. This was confirmed
from secondary sources quoting the construct and buffer, and independently by the reviewer in
[`post_pivot_review.md`](post_pivot_review.md) §2, which reached "pH 6.3 on a
C6A/C111A/F50E/G51E background" from the full text. **The primary table itself was not
accessible to me** — see §6.

## 3. Finding 1 — the disulfide mismatch (the important one)

| | experiment | our simulation |
|---|---|---|
| metal | apo | apo ✅ |
| oligomer | monomer (F50E/G51E) | monomer ✅ |
| **Cys57–Cys146** | **intact (SS)** | **reduced (2SH)** ❌ |
| free cysteines | removed (C6A/C111A) | **present** (wild-type C6, C111) ❌ |
| pH | 6.3 | implicit ~7 (no titratable-state model) ❓ |
| temperature | 25 °C = 298.15 K | 298.15 K ✅ |

The reduced state is not an accident — it is rule 1 of `CLAUDE.md` and is enforced twice, by
`fep.keep_disulfide_reduced` answering pdb2gmx's `-ss` prompts and by
`assert_topology_disulfide_free`. `README.md` §68 justifies it as the disease-relevant species
and notes Wells 2021 and Hsueh 2022 simulated apo-2SH.

**That justification is about biological relevance. It is not a statement that the controls
were measured in that state — and they were not.** The two questions were never separated.
Apo-2SH is a defensible thing to simulate; it is not the thing Lindberg, Nordlund & Oliveberg
or Stathopulos measured.

### How much does this matter?

Honestly: **unknown, and it should not be overclaimed.** ΔΔG is a difference of differences,
so a disulfide contribution that is identical in wild-type and mutant cancels exactly. For
mutation sites far from the 57/146 bond, much of it plausibly does.

What breaks that argument is that the disulfide's effect on SOD1 is not a constant offset —
reduction is known to loosen the β-barrel and the zinc/electrostatic loops, and the folded-state
ensemble of apo-2SH is the *more* flexible one. A mutation's cost depends on the rigidity of
the state it is made in. So the cancellation is an assumption, not a given, and it is the
folded leg — the one this project has never got to converge — where it would fail.

Two observations bear on it, pointing in opposite directions:

- **Against the disulfide as the main error:** the gate's residuals are not a uniform offset. Non-glycine errors are all positive (+0.46, +0.93, +1.34, +1.74) and all three position-93 errors negative (−1.17, −2.50, −4.17). A reference-state difference common to all eight points cannot produce a site-dependent sign flip. It cannot be the whole story, and it is not an alternative explanation for the position-93 compression.
- **For it mattering:** it is a systematic, un-modelled difference affecting every point, and unlike sampling it **cannot be reduced by spending GPU hours.** Gate attempt 2 already established by direct test that sampling was not the limiting error (A4V, 3 ns → 9 ns folded: precision 5×, accuracy 9%). This is a candidate for what the residual actually is, and it is testable.

**The cheap test exists.** `keep_disulfide_reduced` is a config flag, and the engine already
has a guard that asserts the topology is disulfide-free. Running one gate variant in the
oxidized state on the SS-matched reference would measure the cancellation directly rather than
assuming it. That is a scope change (rule 1 forbids switching states without sign-off), it is
one variant of GPU time, and it is the highest-information experiment currently available to
this project. **Not started — it needs the user's decision.**

## 4. Finding 2 — the Lindberg-derived values are exactly right

Lindberg 2005 Table 1 reports ΔG per mutant on the monomeric apoSOD background, not ΔΔG.
Computing ΔΔG = ΔG(pWT) − ΔG(mutant):

| variant | ΔG (Table 1) | ΔΔG computed | `variants.csv` | |
|---|---|---|---|---|
| pWT | 3.03 ± 0.11 | — | — | reference |
| A4V | 1.41 ± 0.05 | **1.62** | 1.62 | ✅ exact |
| G93A | 0.60 ± 0.05 | **2.43** | 2.43 | ✅ exact |
| I113T | 1.78 ± 0.04 | **1.25** | 1.25 | ✅ exact |

Three independent values reproducing to the stated precision is strong evidence the panel's
monomer column was derived correctly and from the right table, with the right sign convention
(positive = destabilizing) and the right background. Whatever is wrong with the gate, **it is
not a transcription error in these three rows.**

Note this also settles a smaller worry: `exp_ddg` and `exp_ddg_dimer` are not swapped. I113T's
1.25 (monomer) and 2.48 (dimer) are distinct quantities from distinct tables.

## 5. Finding 3 — `measured_state: apo_monomer` is likely wrong for G93S and G93V

All eight panel rows carry `measured_state: apo_monomer`. For the Stathopulos pair this looks
incorrect on two counts:

1. **Oligomer.** Stathopulos's pseudo-WT is described as "a pseudo wild-type background containing no free cysteines" — i.e. C6A/C111A. No F50E/G51E interface mutations are reported. Without them apo SOD1 is **dimeric**, so these are apo-*dimer* measurements, not apo-monomer.
2. **Method and temperature.** DSC measures thermal unfolding; destabilization is reported as ΔTm (≈8 °C for G93S, ≈16 °C for G93V). Converting that to a ΔΔG at 25 °C requires a ΔCp model and extrapolation from the Tm (≈50–60 °C). That is a different observable from Lindberg's urea ΔΔG at 25 °C, obtained by a different route.

This is very likely the actual explanation of a long-standing oddity in this repo — recorded in
[`f64a_folded_leg_failure.md`](f64a_folded_leg_failure.md) as "G93V and G93S … the only
controls where monomer > dimer, all from Stathopulos 2006", and flagged in `HANDOFF.md` §5.1 as
"three experimental values are suspect". **They are not suspect measurements. They are
correctly measured values of a different quantity, filed in the monomer column.**

That matters directly: G93V is the largest single contributor to the gate failure (17.38 of
31.56 SSE), and G93S is second. **The two worst gate points are the two whose reference state
is least comparable.** That is not a reason to drop them — dropping inconvenient controls is
the post-hoc move `CLAUDE.md` forbids — but a gate that mixes urea-monomer and DSC-dimer
references is measuring reference heterogeneity as well as calculation error, and cannot
separate the two.

## 6. Finding 4 — two repo statements about the sources are wrong

| where | says | actual |
|---|---|---|
| `citations.md:11-14`, `:145` | `Kumar2017_TableS1` — "**UNVERIFIED — not found on PubMed**", "confirm this citation before publishing anything that rests on it" | **It exists.** Sci Rep 7:4678, [PMC5498623](https://pmc.ncbi.nlm.nih.gov/articles/PMC5498623/). `README.md:492` had the correct DOI the whole time; the two files disagreed and nobody reconciled them. |
| `config/pipeline.yaml:49` | `controls_csv: … # verified apo-monomer controls (Kumar2017 Table S1)` | "Verified" is wrong in the other direction. The source exists, but the values are a **compilation with no stated conditions**, and at least two rows are not apo-monomer (§5). |

Both are corrected in the same commit as this document. Note the pair is instructive: one file
under-claimed and one over-claimed the *same* citation, and the gate ran on it for a month.

## 7. Conflicts for the user to resolve

Per `CLAUDE.md`, surfaced rather than decided:

1. **Simulated state vs measured state.** v1 simulates apo-2SH (rule 1, non-negotiable, "if holo seems needed, STOP and ask"). Every gate control is apo-SS. Either the simulation moves to the state the controls were measured in, or the gate needs controls measured in apo-2SH, or the mismatch is accepted and stated as a limitation in every result. Rule 1 makes this the user's call. **It is the single most consequential open question in the project** — more so than any remaining sampling question, because sampling was already excluded by direct test.
2. **Free cysteines.** The controls are all C6A/C111A. We simulate wild-type C6 and C111. Same class of issue, smaller magnitude, and it interacts with (1): C6A/C111A exists precisely to stop the free thiols scrambling with the 57–146 bond.
3. **pH 6.3 vs our implicit neutral.** Fixed-charge FEP has no titratable states, so this cannot be modelled directly, only bounded. Worth checking whether any gate site has a titratable neighbour whose protonation would plausibly shift between 6.3 and 7.
4. **The `measured_state` column.** If §5 is right, G93S/G93V should be relabelled `apo_dimer` and the gate subset reconsidered — but changing the gate composition after seeing the results is exactly the post-hoc move the pre-registration forbids. The honest options are to keep them and state the heterogeneity, or to pre-register a corrected subset for a *future* gate. **Do not silently relabel and re-evaluate.**

## 8. What is verified, and what is not

Verified against primary or quoted primary text:

- Kumar 2017 exists, is a compilation, names its five upstream sources, states ~0.3 kcal/mol methodological error, and gives no experimental conditions.
- Lindberg 2005: construct, monomeric variant, apo/EDTA, pH 6.3, 25 °C, urea, **disulfide intact** (direct quotation), and the ΔG values that reproduce A4V/G93A/I113T exactly.
- Stathopulos 2006: DSC, cysteine-free pseudo-WT, apo and holo, G93 series, ΔTm magnitudes.

**Not verified — do not cite these as established:**

- **Nordlund & Oliveberg 2006's primary table.** PNAS returned 403 and I could not locate a PMC copy. The construct and buffer are corroborated by secondary sources and by the reviewer's independent reading, but **the ΔG/ΔΔG values behind F64A (−0.20), I149A (4.05) and I18V (0.37) have not been checked against the source table.** F64A is the value this repo has questioned twice; it remains unchecked at the primary level. Highest-priority remaining item.
- **Stathopulos's numeric ΔΔG** for G93S (3.70) and G93V (7.00), and how a ΔTm was converted into them. JBC returned 403.
- Whether Stathopulos used any interface mutation. Absence of a report is not proof of absence; §5's dimer inference rests on that absence.
- Kumar's Table S1 itself — whether our eight values match his rows, and whether he applied any harmonization.
- The Vassall 2006 and Byström entries, which supply other rows of the 54 but none of our eight.

**A note on method.** Several of these facts were extracted by a summarizing fetch rather than
read directly, and one early extraction misreported Lindberg's Table 1 as ΔΔG (−1.62/−2.93/
−2.48) when it holds ΔG (1.41/0.60/1.78). That error was caught only because the ΔG route
reproduced the panel exactly and the ΔΔG route did not. **Treat any single machine-extracted
number in this document as provisional until read from the source PDF.** The user has library
access through BU; the three papers in §6's "not verified" list should be pulled directly.

## 9. Next

1. **Get the three PDFs** (N&O 2006, Stathopulos 2006, Kumar 2017 Table S1) through the BU library and check the seven unverified numbers. Nothing else in this audit can be closed without them.
2. Decide conflict 7.1. Everything downstream — whether the gate can pass, what the LiveCoMS note claims, whether weeks 2–3 are worth running — depends on it.
3. If 7.1 resolves toward matching the experiment, the one-variant SS test in §3 is the cheapest decisive measurement available.
