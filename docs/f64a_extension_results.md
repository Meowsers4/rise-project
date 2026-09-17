# F64A exact-path extension result — 2026-09-16

## Evidence and execution

The authorized folded-only apo-2SH pilot extended the same 20 lambda windows × three
repeats from 3 to 9 ns of retained sampling. Production array **7589742** finished and
the SCC validator accepted all 60 production NPZs: one protocol/code identity, finite
arrays of shape `(20, 9001)`. Successful first-light job: **7589640**. Earlier first-light
attempts **7585581** and **7587677** were preserved separately after, respectively, the
temporary-TPR filename failure and the restart-boundary integrity diagnosis.

The actual SCC source was
`/projectnb/rise-batteries/bode/rise-project/results/fep/F64A/folded`;
the previously documented SCC `archive_2026-09-11` directory was absent. Before staging,
all 60 NPZs matched the frozen forensic manifest, and all 60 XVG histories reproduced
their NPZ energies exactly. The historical results tree was not used as an MD working
directory. See the [registered design and boundary amendment](../pilots/f64a_sampling_extension_v1/PREREGISTRATION.md).

The SCC report supplied by the user is preserved without its shell prompt in
[`f64a_extension_review/f64a_extension.json`](f64a_extension_review/f64a_extension.json).
Simulation/primary-analysis code commit: `80097dcb2ff90ebfdedfadb6e4d9a866a3d9b55a`;
extension protocol: `bf6841ccb3b9de79`. Config and module hashes are in that report.
This local copy is a user-supplied cluster result, **not** a fresh local MBAR reproduction:
the extended raw NPZs remain on the SCC.

## Result

All entries below are **folded-leg mutation free energies**, in kcal/mol; they are
not folding ΔΔGs or experimental-accuracy measurements. The ± values are repeat SEMs,
not 95% confidence intervals. For the cumulative estimates, the report's larger-of-
propagated-MBAR-error-and-repeat-SEM rule selected the repeat SEM in all three cases.
The final-block summary is calculated from its three reported repeat estimates.

| Interval | Mean ΔG ± repeat SEM | Repeat range |
|---|---:|---:|
| Cumulative 3 ns | 1.704 ± 0.402 | 1.256 |
| Cumulative 6 ns | 1.224 ± 0.188 | 0.609 |
| Cumulative 9 ns | 0.771 ± 0.105 | 0.360 |
| Disjoint final 3 ns (6–9 ns) | 0.587 ± 0.100 | 0.344 |

The cumulative mean shifted by **−0.933298 kcal/mol** from 3 to 9 ns. Its successive
3→6 and 6→9 ns changes were −0.479835 and −0.453463 kcal/mol, respectively. The repeat
range fell **71.38%**. Increased duration therefore changed the estimates substantially
as well as narrowing repeat disagreement; it did not merely refine the original value.

| Repeat | 3 ns | 6 ns | 9 ns | Final 3 ns |
|---|---:|---:|---:|---:|
| r0 | 1.250 | 1.091 | 0.607 | 0.558 |
| r1 | 2.506 | 1.596 | 0.739 | 0.429 |
| r2 | 1.357 | 0.987 | 0.967 | 0.773 |

Repeat r1 was the clearest instability: its 3 ns estimate differed from its later
disjoint block by **2.076474 kcal/mol**, despite 3 ns signed-net hysteresis of only
**0.022705 kcal/mol**. Its cumulative 9 ns net hysteresis was 1.104386 kcal/mol;
sum-absolute and maximum-local discrepancies changed from 1.678103/0.222590 at 3 ns to
2.889829/0.961326 at 9 ns. The local diagnostics did not improve uniformly. The other
repeats also moved; their 3 ns-versus-final-block differences were 0.691310 and
0.583723 kcal/mol. At 9 ns the minimum adjacent overlaps were 0.0204–0.0241.

## Bounded interpretation and next step

**Established:** 3 ns was insufficient to obtain duration-stable folded-leg estimates
for these paths. Extension reduced repeat disagreement. Small early aggregate
hysteresis did not certify a stable estimate.

**Not established:** convergence at 9 ns, adequate conformational coverage, a physical
slow-mode mechanism, agreement with experiment, or a new validation-gate pass. All
cumulative estimates are dependent, and each interval re-estimates equilibration
trimming/correlation thinning. The final block is disjoint from the first 3 and first
6 ns, but is contained within the cumulative 9 ns estimate. These are descriptive
comparisons, not independent hypothesis tests.

The user authorized CPU-only follow-up after reading this result. The next analysis
compares disjoint 0–3, 3–6, and 6–9 ns blocks and the two halves of the final block,
while recording per-window trimming, correlation estimates, retained sample counts,
overlap matrices, and signed local discrepancies. See
[`f64a_block_analysis_design.md`](f64a_block_analysis_design.md). No additional GPU
arm is selected or authorized by this report. The completed
[block results](f64a_block_analysis_results.md) show downward movement in all three
repeats within the final block, alongside changing trimming choices. The user then
authorized the isolated [CPU selection sensitivity](f64a_selection_sensitivity_design.md).
The frozen failed historical gate remains
unchanged. This result and the existing benchmark audit should be reviewed by
**Peter Vlasov**; no mentor message has been sent by the agent.
