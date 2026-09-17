# F64A disjoint-block results — 2026-09-16

The user completed the CPU analysis on SCC and supplied the full
[report](f64a_extension_review/f64a_blocks_v1.json). It records analyzer commit
`0ba895ab0c88dbeb11d086b542ea07b79d00ad2a`, simulation commit
`80097dcb2ff90ebfdedfadb6e4d9a866a3d9b55a`, and protocol `bf6841ccb3b9de79`.
All 15 fits, 300 window records and 285 adjacent records are present. The report
records successful reproduction of the original first/final 3 ns fits. Local review
verified all five source/config/reference hashes and summary/contrast arithmetic;
extended matrices remain on SCC, so this is not independent local MBAR reproduction.

All values are folded-leg mutation free energies, not folding ΔΔGs. Errors below
are repeat SEMs, not confidence intervals or uncertainty covering unsampled modes.

| Disjoint interval | Mean ΔG ± repeat SEM, kcal/mol | Repeat range |
|---|---:|---:|
| 0–3 ns | 1.704 ± 0.402 | 1.256 |
| 3–6 ns | 1.223 ± 0.153 | 0.499 |
| 6–9 ns | 0.587 ± 0.100 | 0.344 |
| 6–7.5 ns | 0.953 ± 0.127 | 0.383 |
| 7.5–9 ns | 0.404 ± 0.072 | 0.218 |

The final two halves moved downward in every repeat:

| Repeat | Earlier half | Later half | Signed movement, kcal/mol |
|---|---:|---:|---:|
| r0 | 0.826 | 0.473 | −0.353 |
| r1 | 0.826 | 0.260 | −0.566 |
| r2 | 1.208 | 0.478 | −0.730 |

Mean paired movement is −0.549472 kcal/mol, paired SEM 0.109218. This is descriptive:
disjoint raw columns do not guarantee independent slow modes. Improved repeat
agreement does not establish stationary late estimates.

Selection sensitivity remains unresolved. Of 120 final-half window fits, 15 discard
over 75% of their interval. Maximum discarded fractions are 92.9% in the earlier
half and 90.6% in the later half. Some fits retain only 33–38 selected observations.
Several windows' cutoff fractions change by more than 50 percentage points between
halves. Correlation/sample-count estimates are heuristic, not independent protein-basin counts.

Pairs 13–14 through 18–19 account for 68.2% and 61.6% of total absolute hysteresis
in the earlier/later halves respectively (pooled across repeats). This localizes
diagnostic disagreement, not the causal source of endpoint movement. No structural
slow mode can be identified from energy records alone.

Conclusion: 9 ns convergence remains unestablished. The next authorized step is the
[CPU selection-policy sensitivity](f64a_selection_sensitivity_design.md), not a new
GPU arm. That sensitivity has now completed; see
[its result](f64a_selection_sensitivity_results.md) and the
[draft 9→15 ns design](f64a_continuation_9_to_15_design.md).
The historical failed gate and F64A exclusion remain frozen. Mentor review is
deferred at the user's direction, not a prerequisite; no mentor message has been sent.
