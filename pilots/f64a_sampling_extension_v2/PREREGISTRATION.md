# F64A 9→15 ns bounded continuation — registered 2026-09-17 before execution

The user selected implementation and pre-submission work on 2026-09-17 after
reviewing the completed v1 sensitivity report. This is exploratory development on
known unstable paths, not a blind prediction, convergence gate or experimental test.
The [design](../../docs/f64a_continuation_9_to_15_design.md) defines the question and
interpretation boundaries; this registration fixes the executable contract.

## Fixed scientific arm

Same F64A folded monomer, apo-2SH. Continue every 20 lambda states × three repeats
from v1's 9500 ps checkpoints to 15500 ps, adding exactly 6 ns per path (360 ns
aggregate) and reaching 15 ns retained sampling. No rebuilding, new protein starts,
new repeats, altered dynamics or regenerated velocities. Extend v1 `extension.tpr`,
not the original 3500 ps TPR. All continuation files are isolated below
`results/pilots/f64a_sampling_extension_v2/`; v1 and historical inputs stay unchanged.
Energy/checkpoint/log retention only; structural mechanisms are not an endpoint.

Before staging, all 60 v1 NPZ hashes must match the preserved sensitivity report.
Validate one protocol/identity, `(20, 9001)` finite matrices, original source ancestry,
the actual TPR end time/lambda/energy cadence/no-coordinate settings and checkpoint
time. XVG must have an exact integer-ps grid. Its retained history must equal v1's
frozen matrix except the 3500 ps column, whose maximum difference must match the
existing v1 NPZ provenance exactly. Compare all cluster reports with frozen references.

Hash every source file; capture clean committed code/config/environment/reference
dependencies, ancestry, and resource approval in manifests. First light is a separate
10 ps w0/r0 continuation from the exact 9 ns source. It must pass before production
staging/submission. Production uses fresh source copies, never first-light outputs.

During append, retain the original recorded 3500 ps exception and permit a new
recorded change only at 9500 ps. Any other historical change is fatal. Keep all frozen
9 ns columns unchanged in output NPZs, admitting only 9501…15500 ps new records.
Production shape is `(20, 15001)`; first-light shape `(20, 9011)`. No gaps, duplicates
or non-finite values. Validate existing outputs fully, not by shape alone. Use atomic
TPR creation and validated NPZ publication. Relevant code/config is frozen during jobs.

## Fixed analysis

`analysis.yaml` fixes five blocks: 6–9, 9–12, 12–15 ns, and the two 1.5 ns final
halves. Each uses adaptive, fixed 0%, and fixed 25% trimming with correlation thinning:
45 fits. Supporting adaptive cumulative 9/12/15 ns estimates add nine fits, for 54
total. Primary endpoint is paired 9–12→12–15 ns movement under each policy. Secondary
endpoints are 6–9→9–12 ns and final-half→final-half movement. No result is omitted.

Report all repeat movements, mean/paired SEM, within-block repeat SEM/range,
conditional MBAR/larger-of uncertainties, complete overlap, net/sum-absolute/max-local
hysteresis, selected counts/cutoffs/g/times, solver notes and adjacent diagnostic
changes. Fixed policies estimate g on the fixed suffix, without adaptive detection.
Reproduce v1 cumulative 9 ns and adaptive 6–9 ns metrics/counts within 1e-8 numerical
tolerance. Hash all inputs before/after. Refuse mixed/incomplete inputs, failed fits,
non-finite results and overwrites. CPU analysis cannot schedule GPU producers.

Policies/cumulative fits share data; disjoint blocks can remain slow-mode correlated.
Sample counts/uncertainties are conditional, not guarantees of basin coverage.
Adjacent changes are not a causal decomposition of endpoint drift. No new convergence
pass threshold, folding ΔΔG, gate pass or experimental accuracy claim is introduced.

## Resource admission and finite stop

One L40S/eight CPU slots per task, at most eight simultaneously; 12 hours is a
per-attempt upper limit. Before staging, `admit` requires measured accounting for
job 7589742, project-quota and balance evidence, available GiB, balance SU, applicable
SU/core-hour rate, explicit aggregate GPU-hour/SU caps, and a finite maximum number
of attempts per task. No balance, quota headroom, node rate or budget is guessed.

Accounting must cover exactly 60 successful v1 tasks. Apply the config's conservative
2× peak-runtime margin. Reserve one hour per first-light attempt; derive a production
walltime no greater than 12 hours such that all allowed task attempts plus first light
fit the approved GPU-hour cap. Require that walltime covers measured peak runtime
with the declared margin. Multiply reserved hours by eight slots and the explicit
SU rate, and require both the approved SU cap and current balance to cover it.
This is a conservative reservation bound, not a forecast of actual expenditure.

Require quota headroom for frozen full staging, growing private append sets, new
matrices and a declared 2× storage margin. The recorded evidence/entered values need
human inspection; filesystem-wide `df` free space is not project quota headroom.
Recheck quota and balance before production if the admission snapshot is stale.

Use only the guarded `submit` CLI, which overrides scheduler walltime with the
admitted cap and records job IDs. Each task checks recorded job identity and uses a
locked finite attempt ledger. Only known unusable-GPU placement failures return 99
for rescheduling, and only while attempt budget remains. Other failures stop loudly.
Technical checkpoint resumes consume remaining attempts and cannot add sampling.
No duplicate array submission, unlimited retry loop or deletion of failure evidence.

**Stop at 15 ns and report the declared analysis regardless of outcome.** Any further
sampling or different-start experiment is a separate design/user decision. Mentor
review is deferred at user direction. Historical gate/exclusion/state remain frozen.

TPR/checkpoint inspection uses the officially documented
[`gmx dump` interface](https://manual.gromacs.org/documentation/2025.3/onlinehelp/gmx-dump.html).
Duration-only edits use [`gmx convert-tpr -extend`](https://manual.gromacs.org/documentation/2025.3/onlinehelp/gmx-convert-tpr.html),
with dumped model/state equality checked after excluding only `nsteps`.
