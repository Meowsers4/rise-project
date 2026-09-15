# Methods-and-limitations deliverable — evidence-locked outline

Status: working assembly plan, 2026-09-14. This is the post-pivot deliverable; it is not a
proposal for another validation attempt.

## Central claim

A pre-registered equilibrium-window pmx/AMBER validation gate on apo-2SH SOD1 failed, and the
post-mortem shows why routine convergence and benchmark hygiene checks were insufficient:
within-ladder hysteresis did not order independent-box disagreement, additional sampling
improved precision without materially improving accuracy in the tested variant, and a
compiled experimental column flattened incompatible redox, oligomeric, construct, and
temperature states. One pre-registered SS diagnostic was a converged negative result, so the
redox mismatch is real but does not explain G93A's error.

## Results sequence

### Result 1 — pre-registered gate failure

- Seven usable points; F64A excluded only by the pre-registered closure rule.
- Pearson r = 0.326; RMSE = 2.123 kcal/mol; median closure = 0.220/0.225 depending on rounded
  versus full-precision presentation.
- Preserve the wording “gate of record” and the pivot trigger. Do not re-score it.
- Evidence: [`gate_attempt_1.md`](gate_attempt_1.md),
  [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md).

### Result 2 — the standard diagnostic did not identify the bad replicate

- 48 (variant, leg, replicate) records.
- Folded hysteresis versus independent-box disagreement: Pearson +0.072, Spearman +0.080,
  n = 24.
- F64A folded r1: hysteresis 0.0227 (rank 10/48) and disagreement 1.2026 kcal/mol (rank 1/48).
- Claim only non-separation in this archive; do not generalize to zero information in all FEP.
- Evidence: [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) §6.

### Result 3 — more sampling fixed precision, not accuracy

- A4V folded sampling 3→9 ns/window; unfolded leg reused unchanged.
- Folded box spread 1.46→0.27 kcal/mol; ΔΔG error +1.92→+1.74 kcal/mol.
- Report “precision improved about fivefold, accuracy by 9%.”
- This excludes additional sampling as the useful next action for A4V, not for every variant.
- Evidence: [`gate_attempt_2_abandoned.md`](gate_attempt_2_abandoned.md).

### Result 4 — the experimental benchmark was heterogeneous

- All eight values reproduce primary tables.
- All are apo-SS experiments; calculations are apo-2SH.
- Six are 25 °C monomeric urea measurements on a dimer-splitting construct.
- G93S/G93V are 49.4 °C whole-dimer DSC values on C6A/C111S, mislabeled apo-monomer by the
  compilation.
- Evidence: [`reference_state_audit.md`](reference_state_audit.md).

### Result 5 — reference choices dominate the aggregate metric

- Frozen seven-point predictions and frozen F64A exclusion.
- Halving G93S/G93V alone: r = 0.657, RMSE = 1.122.
- 25 °C correction alone: r = 0.539–0.642, RMSE = 1.187–1.401 across ΔCp models.
- Both: r = 0.754–0.767, RMSE = 1.074–1.117.
- These are not gates. The 25 °C source uncertainty is about ±2.4 kcal/mol per whole-dimer
  ΔΔG, and halving does not create a monomer experiment.
- Evidence and executable figure:
  [`reference_sensitivity_analysis.md`](reference_sensitivity_analysis.md),
  [`reference_sensitivity.py`](../src/analysis/reference_sensitivity.py).

### Result 6 — a pre-registered negative disulfide test

- G93A folded ΔG shift, SS−2SH = −0.0169 kcal/mol.
- Pre-declared negligible band: |Δ| < 0.3 kcal/mol.
- 120 SS windows complete; all three folded topologies physically contain C57–C146.
- Interpretation: disulfide mismatch does not explain G93A's error; no claim about all sites.
- Evidence: [`prereg_g93a_disulfide_diagnostic.md`](prereg_g93a_disulfide_diagnostic.md).

## Figure and table package

| item | content | source/status |
|---|---|---|
| Figure 1 | gate-of-record predicted vs experimental scatter, F64A exclusion identified | `figures/gate_of_record.png`; executable package complete |
| Figure 2 | folded hysteresis vs independent-box disagreement, all 24 folded records | `figures/convergence_diagnostic.png`; executable package complete |
| Figure 3 | A4V 3 ns vs 9 ns replicate/accuracy comparison | `figures/a4v_sampling_sensitivity.png`; executable package complete |
| Figure 4 | four-panel reference sensitivity | `docs/figures/reference_sensitivity.png`; executable now |
| Figure 5 | G93A 2SH vs SS folded ΔG by replicate, topology inset/table | `figures/g93a_disulfide_diagnostic.png`; executable package complete |
| Table 1 | protocol history and hashes, per leg | `raw_result_reconciliation.md` |
| Table 2 | eight-control construct/state/source audit | `reference_state_audit.md` |
| Table 3 | failed gate plus all forensic sensitivity rows | `reference_sensitivity_analysis.md` |
| Table 4 | limitations mapped to what was and was not tested | draft below |

## Limitation table

| limitation | what the evidence supports | what it does not support |
|---|---|---|
| one target / one force field | reproducible SOD1 failure record | general force-field ranking |
| seven usable mutations at five sites | localizes leverage to G93 references | population-level accuracy |
| A4V sampling test only | additional sampling did not close A4V error | sampling excluded for every mutation |
| G93A SS test only | SS does not change G93A mutation cost materially | redox cancellation is universal |
| dimer-derived sensitivity | benchmark normalization materially moves r/RMSE | halved values are monomer measurements |
| no retained trajectories | energy-space diagnostics are reproducible | structural mechanism for the residual |
| A4V attempt-1 windows overwritten | run-time record preserves the point | full raw reanalysis of all seven points |
| no charge-changing implementation | no claim on those variants | C1 coverage or charge-changing VUS triage |

## Exact methods facts to preserve

- GROMACS + pmx, `amber99sb-star-ildn-mut`, equilibrium windows, MBAR.
- 20 lambda states, three independently solvated replicates, folded and capped-tripeptide
  unfolded legs.
- Attempt-1 production/equilibration: 3.0/0.5 ns per window for both legs.
- Protocol `822108e9db71124d` for the attempt-1 archive.
- Thresholds fixed 2026-08-07: r ≥ 0.70, RMSE ≤ 1.5 kcal/mol, median closure ≤ 0.75,
  minimum six usable points; per-variant closure cap 1.0.
- A4V original point is attested, not raw-window-rederived; state this wherever full
  reproducibility is claimed.
- The SS and 2SH G93A arms share a hash; topology + manifest carry redox provenance.

## Forbidden overclaims

- No “first,” “novel FEP protocol,” or “we show SOD1 variants are destabilizing.”
- No “cycle closure is useless”; say it did not separate independent-box disagreement here.
- No mechanistic attribution of the G93 compression to glycine entropy without trajectories.
- No claim that the experiments are wrong.
- No claim that the combined sensitivity row passes validation.
- No claim that SS mismatch explains the gate failure.

## Remaining closeout work

1. Have a second person reproduce Table 3 and Figure 4 from the archived inputs.
2. If a submission-ready artifact is wanted, convert this evidence-locked outline into a
   polished report or manuscript. The LiveCoMS inquiry was parked by user direction on
   2026-09-14 and is not a prerequisite for closing the technical project.

SCC archival verification was completed on 2026-09-14: the live G93A tree is the 120-window
2SH baseline, and the SS diagnostic has an identical second copy with its manifest at
`/projectnb/rise-batteries/bode/archive_2026-09-14/G93A_SS_diagnostic`.

Figures 1--3 and 5 are built by [`methods_figures.py`](../src/analysis/methods_figures.py);
the archived-input contract and exact reproduction command are in
[`methods_figure_package.md`](methods_figure_package.md).
The raw and derived inputs are covered without duplication by
[`forensic_archive_manifest.json`](../data/forensic_archive_manifest.json).
