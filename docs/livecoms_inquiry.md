# LiveCoMS presubmission inquiry — draft, 2026-09-12

Week-1 deliverable from [`post_pivot_review.md`](post_pivot_review.md) §3 ("requires a
presubmission inquiry — **send it in week 1**"). **Draft for the user to review, edit and
send** — not sent, and not to be sent without checking the two items in §3 below.

Target: *Living Journal of Computational Molecular Science*, **Lessons Learned** article type,
which explicitly accepts failed studies and reproducibility investigations.

---

## 1. The draft

> **Subject:** Presubmission inquiry — Lessons Learned: a pre-registered alchemical FEP
> validation gate that failed, and what the standard diagnostics did not predict
>
> Dear Editors,
>
> I would like to ask whether the following would be in scope for a *Lessons Learned* article.
>
> I ran a pre-registered validation gate for alchemical free-energy calculations (GROMACS +
> pmx, AMBER99SB\*-ILDN, MBAR) on eight ALS-associated SOD1 variants with published
> experimental stabilities, in the apo, disulfide-reduced monomer. The acceptance thresholds
> — Pearson r ≥ 0.70, RMSE ≤ 1.5 kcal/mol, median |cycle closure| ≤ 0.75 kcal/mol — were fixed
> before any variant was evaluated. The gate failed: r = 0.326, RMSE 2.12 across seven usable
> points.
>
> What I think is worth reporting is not the failure but three things the post-mortem
> established, each of which cost real compute or real reading to learn:
>
> **1. The convergence diagnostic did not carry the information it is usually trusted for.**
> Across 48 (variant, leg, replicate) records, within-ladder forward/reverse hysteresis is
> near-uncorrelated with whether a replicate agrees with independently solvated repeats of
> itself (folded leg: Pearson +0.07, Spearman +0.08, n = 24). The single sharpest case is a
> replicate whose hysteresis sits in the lowest quintile of the whole dataset while carrying
> the largest replicate disagreement in it (1.20 kcal/mol). I want to be careful here: this
> does not identify which replicate is correct, and two agreeing replicates are not ground
> truth. The claim is only that the diagnostic did not separate them.
>
> **2. Sampling was excluded as the limiting error by direct test, not by argument.** Tripling
> the folded-leg sampling (3 → 9 ns/window) on one variant improved precision about fivefold —
> box-to-box spread 1.46 → 0.27 kcal/mol — and accuracy by 9%. That is a negative result that
> saved roughly 300 GPU-hours, and it is the kind of thing that is rarely published and
> therefore repeatedly rediscovered.
>
> **3. The largest identified problem was a reference-state mismatch, found by reading the
> primary sources rather than by computing.** Every experimental control had been taken through
> a compilation paper that states no experimental conditions. Going back to the primary
> measurements, all of them were made with the Cys57–Cys146 disulfide intact, while the
> calculations were run on the reduced form; and two of the eight were calorimetric
> measurements on the apo dimer, filed in a monomer column. Both were traceable to a single
> compilation step where the conditions were dropped. I suspect this failure mode — an
> experimental benchmark column assembled from a review, with heterogeneous constructs and
> redox states flattened into one number — is not specific to SOD1.
>
> All estimates in the study have been re-derived from the archived per-window reduced
> potentials in three independent local environments spanning both pymbar backends, and agree
> with the original run-time records to ≤ 3e-12 kcal/mol, so the numbers reported are
> reproducible rather than transcribed.
>
> I am aware of the limits. This is one target, one force field, one estimator, one water
> model, and 48 records from eight variants at five sites — not a multi-system study, which I
> understand *Lessons Learned* prefers. The convergence-versus-correctness principle itself is
> prior art; the contribution I am claiming is the reproducible failure record, the
> reference-state audit as a procedure, and the operational lessons, not the principle.
>
> Would this be worth developing into a full submission, or is the single-system scope
> disqualifying? I would rather ask now than write it up against the wrong target.
>
> Thank you for your time.
>
> [name, affiliation, contact]

## 2. What this deliberately does not claim

Checked against `CLAUDE.md`'s forbidden framings and the review's list of withdrawn claims:

- No "first", no "novel FEP protocol", no "we show SOD1 variants are destabilizing".
- Does **not** say the diagnostic contains no information about accuracy — seven points cannot establish that, and the letter says only that it did not separate the cases here.
- Does **not** claim independent solvent boxes are independent conformational sampling.
- Does **not** attribute the position-93 discrepancy to backbone entropy, or to anything.
- Does **not** claim the experimental values are wrong. §4.1 of the audit closed that off: F64A's ΔG is measurably *above* its pWT reference.
- Does **not** present the reference-state mismatch as the explanation of the gate failure. It cannot be — the residuals are site-dependent and a common reference offset is not.

## 3. Before sending — two things I could not settle

1. **Check the current editors and the submission route.** I did not look up who to address or whether LiveCoMS wants inquiries by email or through their GitHub/website. Do not send to a name I did not verify.
2. **The apo-2SH vs apo-SS decision is still open** ([`reference_state_audit.md`](reference_state_audit.md) §7.1). The letter above is written to be true either way, but if you decide to re-run against SS-matched references, point 3 becomes "here is the mismatch and here is what happened when we fixed it", which is a materially stronger letter. **It may be worth deciding that first** — the inquiry is week-1 in the plan, but a two-week delay that converts a limitation into a controlled comparison is probably a good trade.

## 4. Status

Not sent. Awaiting the user's edit and the §3.2 decision.
