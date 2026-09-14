# LiveCoMS presubmission inquiry — ready for user review, 2026-09-14

Week-1 deliverable from [`post_pivot_review.md`](post_pivot_review.md) §3. This is a draft for
the user to complete and send; it has not been sent and must not be sent autonomously.

The current route was verified 2026-09-14 against the official
[author instructions](https://livecomsjournal.github.io/authors/policies/) and
[editorial team](https://livecomsjournal.org/index.php/livecoms/about/editorialTeam): send a
presubmission letter of no more than one page to the Lessons Learned Lead Editor, Daniel M.
Zuckerman, at `lessonslearned@livecomsjournal.org`. A full manuscript is prepared only after
approval; eventual submission is through the journal site.

## 1. One-page draft

> **To:** lessonslearned@livecomsjournal.org<br>
> **Subject:** Presubmission inquiry — Lessons Learned from a failed pre-registered SOD1 FEP
> validation gate
>
> Dear Dr. Zuckerman,
>
> I am writing to ask whether the following would be suitable for a LiveCoMS *Lessons Learned*
> article. We ran a pre-registered validation gate for alchemical free-energy calculations
> (GROMACS + pmx, AMBER99SB\*-ILDN, MBAR) on eight ALS-associated SOD1 variants. The gate
> failed across seven usable points (Pearson r = 0.326; RMSE = 2.12 kcal/mol), and we stopped
> predictive use of the pipeline.
>
> The proposed article would document three linked lessons. First, across 24 folded-leg
> records, within-ladder hysteresis was nearly uncorrelated with independent-box disagreement
> (Pearson +0.07; Spearman +0.08). The worst disagreement (1.20 kcal/mol) occurred at
> hysteresis in the lowest ~20% of all 48 records. Our claim would be system-specific, not that
> hysteresis is generally useless.
>
> Second, tripling folded-leg sampling from 3 to 9 ns/window for A4V reduced box-to-box spread
> from 1.46 to 0.27 kcal/mol but reduced its error by only 0.18 kcal/mol. That negative result
> prevented six further arrays: precision improved about fivefold while nearly all of the
> accuracy error remained.
>
> Third, the benchmark was heterogeneous. All eight experiments retained the Cys57-Cys146
> disulfide while the calculations used apo, reduced SOD1. Six were monomeric urea measurements
> at 25 °C; two values compiled as apo-monomer data were whole-dimer DSC measurements at
> 49.4 °C. A pre-registered SS-versus-2SH G93A test was negative (folded shift -0.017 kcal/mol),
> so redox did not explain that error. The lesson is benchmark chain-of-custody, not that the
> primary experiments are wrong.
>
> Seven estimates re-derive from archived reduced potentials; the eighth survives only as its
> run-time record because resubmission overwrote its windows. We will provide executable
> analyses, checksums, protocol hashes, uncertainty and topology checks, and the failed gate
> unchanged. The contribution is the reproducible failure record, reference audit, and
> compute-stopping decision—not a novel FEP method. It is not adapted from a previous article.
>
> The study covers one target, one force field, one water model, and eight variants at six
> sites. We think its operational lessons apply to simulation studies using small validation
> panels and compiled experiments. Does this meet the category's broad-usefulness requirement,
> or would you recommend another venue or a broader dataset?
>
> The proposed author(s), [names], have expertise in [edit to be accurate]. We would maintain a
> public GitHub repository, review it annually and after material updates, respond to community
> issues, and transfer maintenance to [named coauthor, lab, or organization] if needed. We
> propose a CC BY 4.0 license.
>
> Thank you for your consideration,
>
> [name, affiliation, ORCID, contact]

## 2. Required elements and claim controls

The draft now includes every element listed in the current presubmission instructions: scope,
difference from existing work, prior-article status, author expertise, update and stewardship
plan, and proposed license. The user must replace the bracketed expertise, stewardship, and
identity fields before sending.

It deliberately does not claim:

- “first,” “novel FEP protocol,” or “we show SOD1 variants are destabilizing”;
- that hysteresis contains no information in general, or that agreeing boxes are ground truth;
- that additional sampling is excluded for variants other than the tested A4V case;
- that the primary experiments are wrong or that the disulfide mismatch explains the gate;
- that halving whole-dimer measurements creates isolated-monomer measurements;
- that the post-verdict sensitivity scenarios re-open or pass the gate.

## 3. Verified route and remaining user decisions

Verified 2026-09-14:

1. *Lessons Learned* explicitly includes negative-result and reproducibility articles, but asks
   for value to a substantial subset of the community and prefers multiple systems/conditions.
2. A presubmission letter is mandatory, must be no more than one page, and normally receives a
   response within two weeks.
3. The current Lead Editor is Daniel M. Zuckerman and the route is
   `lessonslearned@livecomsjournal.org`.
4. After approval, authors prepare the article in the LiveCoMS LaTeX template in a public
   GitHub repository and later submit through the journal website. Do not build or submit that
   manuscript before the inquiry is approved.

Before sending, the user still needs to:

1. replace all bracketed identity, expertise, and maintenance-stewardship fields;
2. confirm the bounded scope statement: completed v1 remains apo-2SH with no more GPU work;
3. decide whether the proposed annual/update-triggered maintenance plan and CC BY 4.0 license
   are acceptable;
4. trim or reformat after inserting identity details if the rendered letter exceeds one page.

## 4. Status

**Not sent.** The external route is verified; only the explicit user fields and decisions above
remain.
