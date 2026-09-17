# Draft F64A 9→15 ns exact-path continuation design — 2026-09-16

## Status and question

The user requested documentation and an explanation of the next design. This is
a **proposal**, not an executable package or GPU submission authorization. No new
job has been launched. Mentor review is deferred at the user's direction and is
not a prerequisite. Implementation/submission require explicit user selection;
commit the final registered design and package before execution.

Question: **does another fixed 6 ns of sampling reduce late time dependence in the
same F64A folded paths, across the three declared selection policies?**

The [completed sensitivity](f64a_selection_sensitivity_results.md) found mean
late movement of −0.485 to −0.549 kcal/mol across policies, robust downward
movement in r1/r2, and policy-sensitive r0 movement. Six additional ns permits
two new equal-length 3 ns blocks comparable with the preceding block. It is a
practical bounded choice, not an inferred relaxation time or prediction that
15 ns will converge.

## Simulation arm

- Same folded monomer, apo, disulfide-reduced (apo-2SH), F64A mutation.
- All 20 lambda states × three repeats from completed v1 production outputs,
  protocol `bf6841ccb3b9de79`: exact 9 ns checkpoints, TPRs, append files and
  frozen energy matrices. No rebuilding, new protein starts, regenerated
  velocities or new repeats.
  The continuation TPR is v1's production `extension.tpr`, ending at 9500 ps,
  not the original 3500 ps `prod.tpr`. Verify its end time before adding 6000 ps.
- Extend every path by 6000 ps, from simulation time 9500 to 15500 ps;
  retained sampling grows from 9 to 15 ns. Additional aggregate sampling: 360 ns.
- Preserve timestep, temperature, lambda schedule, force field, thermostat/
  barostat, checkpoint state and energy stride. Only TPR duration changes.
- One L40S and eight CPU slots per task, at most eight tasks running
  simultaneously, 12-hour per-attempt walltime, project `rise-batteries`.
- Energy estimates, checkpoints and logs only, matching the original source
  TPR. No structural-mechanism endpoint or new coordinate-output claim.
- Proposed package `pilots/f64a_sampling_extension_v2/`, output root
  `results/pilots/f64a_sampling_extension_v2/`. Never modify the v1 package,
  historical `results/fep/F64A` outputs or base pipeline config to implement v2.

The full ladder remains included because local diagnostics do not establish which
windows cause endpoint movement. This is targeted to one variant and leg, not a
campaign-wide rerun or a different-starting-conformation experiment.

## Integrity and first light

Before staging, match every source NPZ to the hashes in the preserved v1 sensitivity
report. Verify exact shape `(20, 9001)`, finite energies, state/window/repeat metadata,
one simulation identity/protocol, source checkpoint times and XVG histories.
The source XVG carries the documented v1 restart exception at 3500 ps: validate
that exception against existing recorded provenance and require agreement with
the frozen v1 matrix at every other retained time. No new unexplained historical
exceptions are admitted.

Copy and hash **all** source files, not just NPZs; record source ancestry plus the
new package Git/config/module hashes in a stage manifest and every output. GROMACS
must pass checkpoint append checksums. Never run MD in the source directory.

First light: isolated 10 ps continuation of w0/r0 from the 9500 ps checkpoint.
It must pass GPU execution, append compatibility, exact new time grid, energy
parsing, history reconciliation and metadata checks before production. Production
starts from fresh copies of the frozen 9 ns sources, not first-light outputs.

Preserve the v1 energy matrix intact. After continuation, a newly recorded restart-
boundary change may be accepted **only at 9500 ps**, carrying forward the existing
documented 3500 ps exception. Every other historical energy must remain exact.
Admit new columns only at 9501, 9502, …, 15500 ps, without gaps, duplicates or
non-finite values. Combined arrays must be `(20, 15001)` with their first 9001
columns exactly equal to the frozen source. Never replace the source endpoint
with the regenerated 9500 ps row.

Implement/test this boundary policy before first light; do not widen it after a
production mismatch. Preserve failures for diagnosis and require first light to
pass before the array. Freeze relevant code/config paths while tasks are in flight.
Missing/failed tasks prevent complete analysis; never drop them. Technical resume
uses checkpoints without extending beyond the registered endpoint.

## Analysis fixed before execution

Retained raw column i corresponds to simulation time 500+i ps. Half-open slices:

| Block | Column slice | Count/window |
|---|---|---:|
| Prior 6–9 ns | `[6001:9001]` | 3000 |
| New 9–12 ns | `[9001:12001]` | 3000 |
| New 12–15 ns | `[12001:15001]` | 3000 |
| First half of final block | `[12001:13501]` | 1500 |
| Second half of final block | `[13501:15001]` | 1500 |

Solve every block under the three established policies: adaptive trim-and-thin;
zero additional trimming with correlation thinning; fixed first-25% trimming with
correlation thinning. Fixed policies estimate both correlation times on the suffix,
never invoking adaptive equilibration detection. Preserve the historical neighbor
choice and use the larger sampled-state/neighbor-difference g. Do not choose a
policy according to which result appears stable.

**Primary comparison:** paired 9–12→12–15 ns folded-leg movement, reporting each
repeat, repeat mean and paired SEM separately under every policy.

**Secondary comparisons:** 6–9→9–12 ns movement and final-half→final-half movement.
Report all comparisons regardless of sign/magnitude, with within-block repeat
range/SEM, conditional MBAR error and larger-of uncertainty. These five blocks ×
three policies × three repeats require 45 MBAR fits.

Supporting summaries: cumulative 9/12/15 ns under the adaptive policy (nine more
fits), explicitly overlapping/dependent. Total: 54 planned fits. Record aggregate
net, sum-absolute and maximum-local hysteresis, minimum adjacent overlap, full
overlap matrix, solver notes, and all window cutoffs/g/counts/selected times.
Record forward/reverse adjacent diagnostic estimates and paired changes, without
treating them as a causal decomposition of endpoint drift.

Reproduce preserved v1 cumulative 9 ns and adaptive 6–9 ns metrics/counts within
the existing 1e-8 numerical tolerance, not a convergence threshold. Require unchanged
frozen source columns/input hashes. Record new simulation identity separately from
analysis identity and hash inputs before/after analysis. Refuse mixed/incomplete
protocols and overwrites; no plausible partial report after failed fits. Represent
staging, first light, continuation and analysis as isolated Snakemake rules. CPU
analysis must consume completed outputs without scheduling missing GPU producers.

## Resource check and stop contract

The historical estimate for a 60-path 6 ns extension was approximately 36 GPU-hours,
or 4.5 hours occupied time at concurrency eight, excluding queueing. This is **not**
measured v1 accounting. Before selecting a production budget, inspect `qacct -j
7589742` for actual per-task runtimes/CPU charging, confirm current project balance
and storage headroom, and validate staging/output size estimates. Twelve hours per
task is a safety limit, not expected runtime; one full array's per-attempt reservation
ceiling is 720 GPU-hours. Agree aggregate charging and retry limits before submission.
No current storage or queue capacity is assumed.

**Stop at 15 ns and report the declared analysis.** No automatic extension until a
preferred result appears. No new numerical convergence pass threshold is introduced.
Persistent movement means duration stability remains unresolved. Less movement
across policies/repeats supports improved observed stability, not complete basin
coverage. Stable blocks with repeat disagreement may motivate a separately authorized
different-start experiment. Stable, agreeing blocks cannot distinguish shared
trapping from model/reference bias.

This arm cannot establish a validation-gate pass, experimental accuracy, folding
ΔΔG or a structural slow-mode mechanism. The historical gate and exclusion stay
fixed. Conditional uncertainty cannot cover unvisited modes; this limitation and
the fixed-endpoint approach follow the
[alchemical best-practices guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC8388617/).
