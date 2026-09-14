# Reference-state audit: are the experimental controls comparable to what we simulate? — 2026-09-12, updated 2026-09-14

The last week-1 deliverable of [`post_pivot_review.md`](post_pivot_review.md) §5, and the
question [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) §7 explicitly could not reach:
reanalysis proved the arithmetic is right; it cannot prove the calculation and the experiment
describe the same thermodynamic state.

**Headline: they do not.** Every gate control was measured with the **Cys57–Cys146 disulfide
intact**. We simulate the **disulfide-reduced (2SH)** form, deliberately and by rule. The gate
has been comparing apo-**SS** experiment against apo-**2SH** calculation on all eight points.

A second finding runs the other way: **all eight gate values now trace to a primary table** —
A4V, G93A, I113T from Lindberg 2005; F64A, I18V, I149A from Nordlund & Oliveberg 2006;
and G93S/G93V from Stathopulos 2006. The first six re-derive exactly in their intended
monomer reference. The last two also match their primary table exactly, but that check exposed
a category error in the compilation: 3.70 and 7.00 kcal/mol are whole-dimer apo DSC values at
49.4 °C, not apo-monomer values at 25 °C. The problem is the state and normalization, not a
transcription error. See §5.

This is a **scope conflict for the user to resolve**, per `CLAUDE.md` ("surface the conflict
instead of silently choosing"). No protocol parameter, panel value, or gate verdict has been
changed on the strength of it. See §7.

---

## 1. The chain of custody

`data/variants.csv` cites two sources per gate variant: `Kumar2017_TableS1` plus a primary
paper. **Kumar 2017 is a compilation, not a measurement.**

| | |
|---|---|
| Paper | Kumar V, Rahman S, Choudhry H, Zamzami MA, Jamal MS, Islam A, Ahmad F, Hassan MI (2017). *Computing disease-linked SOD1 mutations: deciphering protein stability and patient-phenotype relations.* Sci Rep 7:4678. |
| DOI / PMC | [10.1038/s41598-017-04950-9](https://doi.org/10.1038/s41598-017-04950-9) · [PMC5498623](https://pmc.ncbi.nlm.nih.gov/articles/PMC5498623/) |
| What it is | A benchmark of 8 stability predictors (PoPMuSiC, FoldX, mCSM, …) against a **compiled** experimental set that Kumar labels **54 apo-monomer** and **33 holo-dimer** ΔΔG values. |
| Where its numbers came from | Quoted: *"These data have been taken from Vassall et al. Nordlund and Oliveberg, Lindberg et al. Stathopulos et al., and Bystrom et al."* |
| Experimental conditions stated | **None.** No pH, temperature, method, construct or redox state for the compiled values. |
| Stated uncertainty | *"the methodological error in the experimental ΔΔG is on the order of up to ~0.3 kcal/mol."* |

Two consequences.

**The uniform `exp_ddg_err: 0.3` on every panel row is Kumar's blanket figure**, not a
per-measurement uncertainty. It is one number covering five labs, two denaturants, two
temperatures and both calorimetric and chemical methods. It should not be read as the error on
any individual control, and a gate RMSE of 1.5 is not "3–5σ" in any meaningful sense.

**The compilation is where the conditions were lost.** Kumar's supplementary Table S1 labels
its columns `ΔΔG(m)` and `ΔΔG(d)`, then defines them as apo-monomer and holo-dimer. Our panel
inherited that labeling. For the Stathopulos rows the first column is not monomeric: both
columns came from a native-dimer calorimetric model (§5).

## 2. Per-source audit

| gate variants | primary source | construct | oligomer | metal | **disulfide** | conditions | method |
|---|---|---|---|---|---|---|---|
| A4V, G93A, I113T | Lindberg 2005 PNAS 102:9754 | C6A/C111A **+ F50E/G51E** | **monomer** | apo (10 mM EDTA) | **intact / oxidized** | 10 mM MES **pH 6.3**, **25 °C** | urea equilibrium + kinetics |
| F64A, I149A, I18V | Nordlund & Oliveberg 2006 PNAS 103:10218 ([PMC1502438](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502438/)) | C6A/C111A **+ F50E/G51E** | **monomer** | apo (10 mM EDTA) | **intact / oxidized** | 10 mM MES **pH 6.3**, **25 °C** | urea |
| G93S, G93V | Stathopulos 2006 JBC 281:6184 ([PMID 16407238](https://pubmed.ncbi.nlm.nih.gov/16407238/)) | pseudo-WT **C6A/C111S**; no F50E/G51E | **dimer** | apo and holo | **intact / nonreducing** | 20 mM HEPES **pH 7.8**; DSC scan | differential scanning calorimetry, dimer N₂ ↔ 2U model |

### Verified quotations

Lindberg 2005, on the redox state — the decisive sentence:

> "All experiments were done under oxidizing conditions, with the intramolecular disulphide
> linkage between C57 and C146 kept intact."

and on construct and buffer:

> "the double mutant C6A/C111A (SODpwt) as pseudo-WT for further mutant analysis" …
> "the monomeric variant F50E/G51E/C6A/C111A" … "the standard buffer was 10 mM Mes (pH 6.3)
> with 10 mM EDTA to maintain the proteins metal free"

Nordlund & Oliveberg 2006 — **verified directly 2026-09-12** via [PMC1502438](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502438/), located through the NCBI ID converter after pnas.org returned 403:

> "The mutants were constructed on a background of C6A/C111A and the dimer-splitting
> substitutions F50E/G51E… The standard buffer was 10 mM Mes (Sigma) at pH 6.3 with 10 mM EDTA
> to maintain the proteins metal-free… All measurements were done at 25 °C… **we have chosen to
> analyze SOD molecules with the C57–C146 disulphide linkage intact.**"

The disulfide finding is therefore stated explicitly in **both** primary sources, in their own
words, not inferred. That closes the central question of this audit.

Stathopulos 2006 — **verified directly 2026-09-14** from the full primary text and tables:

- pseudo-WT is defined as **C6A/C111S**, not C6A/C111A;
- DSC used 20 mM HEPES, pH 7.8, under nonreducing conditions where the native disulfide
  remains intact;
- the apo data were analyzed as native dimer to two unfolded monomers, N₂ ↔ 2U, with the
  concentration dependence and fitted molecularity supporting a dimeric cooperative unit;
- Table 2 reports G93S −3.7 and G93V −7.0 kcal/mol at the common apo `tavg` of 49.4 °C;
  Table 4 reports the corresponding holo values −3.1 and −5.4 at 87.2 °C.

The source convention is `ΔG(mutant) − ΔG(pseudo-WT)`, so negative means destabilizing.
Kumar and this repo store the sign-flipped positive magnitudes.

## 3. Finding 1 — the disulfide mismatch (the important one)

| | experiment | our simulation |
|---|---|---|
| metal | apo | apo ✅ |
| oligomer | monomer for 6 controls; dimer for G93S/G93V | monomer; mismatched for 2 controls ❌ |
| **Cys57–Cys146** | **intact (SS)** | **reduced (2SH)** ❌ |
| free cysteines | C6A/C111A for 6; C6A/C111S for 2 | **present** (wild-type C6, C111) ❌ |
| pH | 6.3 for 6; 7.8 for G93S/G93V | prepared at 7.0 ❓ |
| temperature | 25 °C for 6; common 49.4 °C for G93S/G93V | 298.15 K; mismatched for 2 ❌ |

All three primary sources say so outright — Lindberg's "under oxidizing conditions, with the
intramolecular disulphide linkage between C57 and C146 kept intact" and Nordlund & Oliveberg's
"we have chosen to analyze SOD molecules with the C57–C146 disulphide linkage intact". This is
not an inference from construct design; Stathopulos independently states that the native
disulfide remains intact under its nonreducing DSC conditions.

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

- **Against the disulfide as the main error:** the gate's residuals are not a uniform offset. Non-glycine errors are all positive (+0.46, +0.93, +1.34, +1.92 — the gate-of-record A4V value from `gate_attempt_1.md`, not the attempt-2 rerun's +1.74) and all three position-93 errors negative (−1.17, −2.50, −4.17). A reference-state difference common to all eight points cannot produce a site-dependent sign flip. It cannot be the whole story, and it is not an alternative explanation for the position-93 compression.
- **For it mattering:** it is a systematic, un-modelled difference affecting every point, and unlike sampling it **cannot be reduced by spending GPU hours.** Gate attempt 2 already established by direct test that sampling was not the limiting error (A4V, 3 ns → 9 ns folded: precision 5×, accuracy 9%). This is a candidate for what the residual actually is, and it is testable.

**The cheap test was completed 2026-09-14.** The pre-registered G93A SS diagnostic physically
verified C57–C146 in all three independently built folded systems. Its folded-leg endpoint was
9.7493 ± 0.0547 kcal/mol, only −0.0169 kcal/mol from the frozen 2SH mean of 9.7662. That falls
inside the pre-declared negligible band (`|Δ| < 0.3`): the disulfide mismatch is real but does
not explain G93A's ≈1.17 kcal/mol underprediction. This is one mutation, not proof of universal
cancellation, but it removes the disulfide as the leading G93A explanation.

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

**Nordlund & Oliveberg's three verify too**, and that table reports ΔG *and* ΔΔG, so the
convention can be checked rather than assumed:

| variant | ΔG_D-M,H₂O | ΔΔG reported | ΔG(pWT) − ΔG(mut) | `variants.csv` | |
|---|---|---|---|---|---|
| pWT | 2.87 | — | — | — | reference |
| F64A | **3.07** | −0.20 | −0.20 | −0.20 | ✅ |
| I18V | 2.50 | 0.37 | 0.37 | 0.37 | ✅ |
| I149A | **−1.18** | 4.05 | 4.05 | 4.05 | ✅ |

Same convention as Lindberg (ΔΔG = ΔG(pWT) − ΔG(mutant), positive = destabilizing), internally
consistent to the reported precision. **These six monomer values reproduce exactly.** The
panel's experimental column is numerically sound; §3 and §5 are about what those
numbers mean, not whether they were copied correctly.

### 4.1 F64A is not an outlier in a table — it is a measured ΔG above pWT

F64A's ΔG is **3.07 against pWT's 2.87**. The protein is mildly *more* stable than the
reference, and the −0.20 is that difference, not a rounding artefact or a transcription slip.

This retires the argument in
[`f64a_20window_result.md`](f64a_20window_result.md) §"Additional evidence", which extrapolated
an atom-count trend to predict ~8 kcal/mol of destabilization and concluded "the experimental
value is the outlier". The primary measurement sits in a self-consistent series with its own
pWT reference. Whatever explains the ~7 kcal/mol gap against our FEP number, **"the experiment
is probably wrong" is no longer an available explanation**, and per `CLAUDE.md` it was never a
licence to drop the control.

### 4.2 I149A's apo monomer is not folded at zero denaturant

I149A's ΔG is **−1.18 kcal/mol** — negative absolute stability. The apo monomer is
thermodynamically unfolded in the absence of urea, and its 4.05 ΔΔG is obtained by linear
extrapolation from the transition region back to zero denaturant.

Two consequences. First, Kumar's blanket ±0.3 is certainly too tight for this point: an
extrapolated ΔG through a negative value carries the m-value uncertainty as well. Second, and
more awkwardly, **our FEP computes a folded-state ΔΔG for a variant whose experimental folded
state is not populated.** I149A is the **second**-largest gate reference (4.05, behind only G93V's 7.00) and the second most
accurate FEP point (+0.93 error). Neither observation is grounds to remove it, but a gate
anchored partly on an extrapolation from an unfolded protein is measuring something other than
pure calculation error.

## 5. Finding 3 — Kumar's `apo-monomer` label is wrong for the Stathopulos rows

This is now direct rather than inferred. Stathopulos used pseudo-WT C6A/C111S with no
dimer-splitting F50E/G51E substitutions. The paper fits the calorimetry to a native-dimer
unfolding reaction, N₂ ↔ 2U, and explicitly argues that both apo and holo protein remain dimeric
until thermal unfolding. Kumar nevertheless filed the apo column as `ΔΔG(m)` (apo-monomer)
and the holo column as `ΔΔG(d)` (holo-dimer).

| variant | Kumar `ΔΔG(m)` | Stathopulos apo Table 2 | per monomer stated by Stathopulos | Kumar `ΔΔG(d)` | Stathopulos holo Table 4 |
|---|---:|---:|---:|---:|---:|
| G93S | +3.70 | −3.7 per dimer at 49.4 °C | ≈−1.85 | +3.10 | −3.1 per dimer at 87.2 °C |
| G93V | +7.00 | −7.0 per dimer at 49.4 °C | −3.50 | +5.40 | −5.4 per dimer at 87.2 °C |

The sign change is intentional: Stathopulos uses negative for destabilization, while Kumar and
this repo use positive. The factor of two is not: Stathopulos's discussion expresses the apo
effect as approximately −1.8 to −3.5 kcal/mol **per monomer**, exactly half the dimer-level
Table 2 values.

The values are also not a direct conversion of ΔTm alone. Table 2's pseudo-WT, G93S and G93V
melting temperatures at 0.50 mg/ml are 59.0, 50.7 and 43.1 °C (ΔTm −8.3 and −15.9 °C).
The authors fit the DSC data to the dimer model, use a common heat-capacity model
(`ΔCp = 3.30 kcal mol⁻¹ °C⁻¹` for the reported constant-ΔCp calculation), extrapolate each
ΔG to the common apo `tavg = 49.4 °C`, and then calculate
`ΔΔG = ΔG(mutant) − ΔG(pseudo-WT)`. At 25 °C Table 2 gives −3.5 for G93S and −7.0 for
G93V per dimer; the repo instead contains the common-`tavg` values, +3.7 and +7.0.

This resolves the long-standing "monomer > dimer" oddity. **These are not suspect
measurements. They are correctly transcribed measurements of different metal states, with the
apo dimer total filed as a monomer value.**

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
| `config/pipeline.yaml:49` | `controls_csv: … # verified apo-monomer controls (Kumar2017 Table S1)` | "Verified" is wrong in the other direction. The source exists, but the values are a **compilation with no stated conditions**, and the Stathopulos rows are apo/holo dimer DSC values (§5). |

Both are corrected in the same commit as this document. Note the pair is instructive: one file
under-claimed and one over-claimed the *same* citation, and the gate ran on it for a month.

## 7. Conflicts for the user to resolve

Per `CLAUDE.md`, surfaced rather than decided:

1. **Simulated state vs measured state.** v1 simulates apo-2SH (rule 1, non-negotiable, "if holo seems needed, STOP and ask"). Every gate control is apo-SS. Either the simulation moves to the state the controls were measured in, or the gate needs controls measured in apo-2SH, or the mismatch is accepted and stated as a limitation in every result. Rule 1 makes this the user's call. **It is the single most consequential open question in the project** — more so than any remaining sampling question, because sampling was already excluded by direct test.
2. **Free cysteines.** The Lindberg and Nordlund controls are C6A/C111A; the Stathopulos
   controls are C6A/C111S. We simulate wild-type C6 and C111. Same class of issue, smaller
   magnitude, and it interacts with (1): the substitutions exist precisely to stop the free
   thiols scrambling with the 57–146 bond.
3. **pH mismatch.** Lindberg/Nordlund used pH 6.3, Stathopulos used pH 7.8, and the simulation
   was prepared at pH 7.0. Fixed-charge FEP has no titratable-state model, so this cannot be
   modeled directly, only bounded.
4. **The `measured_state` and normalization fields.** Section 5 proves that G93S/G93V are
   apo-dimer totals at 49.4 °C, not apo-monomer values at 25 °C. Correcting metadata is
   necessary, but changing the gate references or composition after seeing the results is
   exactly the post-hoc move the pre-registration forbids. Preserve the gate-of-record; any
   normalized comparison must be labeled a forensic sensitivity analysis or pre-registered as
   a future gate. **Do not silently relabel and re-evaluate.**

## 8. What is verified, and what is not

Verified against primary or quoted primary text:

- Kumar 2017 exists, is a compilation, names its five upstream sources, states ~0.3 kcal/mol methodological error, and gives no experimental conditions.
- Lindberg 2005: construct, monomeric variant, apo/EDTA, pH 6.3, 25 °C, urea, **disulfide intact** (direct quotation), and the ΔG values that reproduce A4V/G93A/I113T exactly.
- **Nordlund & Oliveberg 2006 (closed 2026-09-12):** construct, buffer, pH 6.3, 25 °C, **disulfide intact** (direct quotation), and ΔG/ΔΔG for pWT, F64A, I18V and I149A — all three reproducing `variants.csv` exactly and internally consistent with the reported pWT.
- **Stathopulos 2006 (closed 2026-09-14):** full primary text, methods and Tables 2/4;
  C6A/C111S construct; absence of F50E/G51E; apo/holo native-dimer model; pH 7.8;
  nonreducing, native-disulfide-intact conditions; ΔTm; ΔCp/extrapolation method; and all four
  Kumar values. G93S 3.70/G93V 7.00 are the sign-flipped apo-dimer totals at 49.4 °C;
  3.10/5.40 are sign-flipped holo-dimer totals at 87.2 °C.
- **Kumar supplementary Table S1 (checked 2026-09-14):** visually checked against the
  Europe PMC supplementary PDF. Its G93S/G93V rows are 3.7/3.1 and 7.0/5.4,
  reference 4, and its note defines the columns as apo-monomer and holo-dimer. The first label
  is incompatible with the cited primary experiment.

**Not verified — do not cite these as established:**

- **`variants.csv:42` labels a row `C111A` (WT-relative, Cys→Ala) while N&O measured `A111C`** — restoring the cysteine onto their C6A/C111A pWT. Same substitution, opposite direction, so the ΔΔG sign convention may be inverted for that row. **Not a gate row and pre-existing**, so nothing was changed; recorded because it is the same reference-frame confusion this audit is about, one level down.
- The Vassall 2006 and Byström entries, which supply other rows of the 54 but none of our eight.

**A note on method.** Several of these facts were extracted by a summarizing fetch rather than
read directly, and one early extraction misreported Lindberg's Table 1 as ΔΔG (−1.62/−2.93/
−2.48) when it holds ΔG (1.41/0.60/1.78). That error was caught only because the ΔG route
reproduced the panel exactly and the ΔΔG route did not.

The N&O extraction has a stronger guarantee: that table reports ΔG **and** ΔΔG, and for all
three mutants ΔG(pWT) − ΔG(mutant) equals the reported ΔΔG to the stated precision *and*
matches `variants.csv`. Three independent agreements are hard to produce by misreading. The
Lindberg extraction is corroborated the same way (computed ΔΔG matching the panel on three
rows). The Stathopulos cross-check is now equally direct: each of the four Kumar G93S/G93V
entries matches the corresponding primary Table 2 or Table 4 value after the documented sign
flip. What fails is Kumar's apo-monomer label, not the numeric transcription.

## 9. Next

1. Decide conflict 7.1 using the now-closed eight-point provenance audit and the completed
   G93A SS diagnostic. Preserve the original gate as the gate-of-record.
2. Run a clearly labeled **forensic sensitivity analysis** using the Stathopulos per-monomer
   normalization and 25 °C values. This may quantify reference heterogeneity but must not be
   called a new gate or used to change the preregistered verdict.
3. Audit the remaining non-gate Kumar rows before using them as a future validation panel;
   the Stathopulos category error shows that column labels cannot be trusted without checking
   each primary source.
