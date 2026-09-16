# F64A exact-path sampling extension pilot — registered 2026-09-15

## Question and motivation

The retrospective temporal analysis found that retained folded-leg estimates were still
moving. F64A folded r1 had the largest 25%-versus-final-quarter movement (1.808 kcal/mol)
and the largest full-duration sibling disagreement (1.203 kcal/mol), despite small signed
net hysteresis at both points. This pilot asks whether continuing the same F64A paths from
3 ns to 9 ns changes the folded-leg estimates and their repeat agreement.

This is an exploratory development experiment on a known failure. It is not a new blind
prediction, a validation gate, or an experiment-accuracy test.

## Fixed design

- System/state: the archived F64A apo, disulfide-reduced (apo-2SH), monomer folded leg.
- Inputs: the exact completed checkpoints and 3 ns NPZ energy matrices from protocol
  `822108e9db71124d`; no system is rebuilt and no velocity is regenerated.
- Sampling: extend every one of 20 lambda states × 3 repeats by 6 ns, producing 9 ns of
  retained production per window. Sixty one-GPU tasks run on L40S GPUs, at most eight at
  once, with a 12-hour per-task limit.
- Isolation: archive inputs are copied and hashed before execution. All continuation files
  and combined NPZs live below `results/pilots/f64a_sampling_extension_v1/`; nothing writes
  to `results/fep/F64A` or to the archive.
- Integrity: before staging, every archived XVG must have the exact 0–3500 ps grid and its
  500–3500 ps reduced energies must exactly reproduce the checksum-frozen NPZ. GROMACS then
  enforces its checkpoint append checksums. New output must have the exact 3501–9500 ps grid.
- Code identity: staging requires committed pilot paths and records the Git commit plus the
  SHA-256 of the pilot module; every NPZ carries both and analysis refuses mixed identities.
- Retention: estimates and checkpoint/log files only. The source TPR saved no coordinates,
  and `convert-tpr` cannot add a trajectory output group. This pilot therefore tests energy-
  space time dependence, not structural mechanism.
- First light: one isolated 10 ps continuation (w0/r0) must finish and parse before the
  60-task production array is submitted.

Approximate cost from the prior 3 ns timings is 36 GPU-hours total and roughly 4.5 hours of
elapsed compute at concurrency eight, excluding queue time. This is a planning estimate,
not a resource guarantee.

## Analysis fixed before execution

For each repeat, recompute folded-leg MBAR and the existing energy-space diagnostics at
cumulative 3, 6, and 9 ns. Recompute equilibration/correlation thinning within each selected
interval. Also solve the disjoint final 3 ns block (6–9 ns).

Report continuously, without fitted pass/fail thresholds:

- folded-leg free energy and conditional MBAR uncertainty;
- signed-net, sum-absolute, and maximum-local hysteresis;
- minimum adjacent overlap;
- absolute movement of the 3 ns and 6 ns cumulative estimates relative to the final 3 ns;
- repeat range and repeat SEM at 3, 6, and 9 ns.

The primary evidence is the folded leg. Any folding ΔΔG formed later by subtracting the old
unfolded leg must be labelled a hybrid-duration secondary summary and must not be compared to
experiment as though it were a new matched benchmark.

## Interpretation boundaries

- Continued movement shows that 3 ns was insufficient for these paths; it does not identify
  a physical slow mode or prove that 9 ns is sufficient.
- Narrower repeat spread does not establish convergence because the three paths retain their
  original, closely related starting ensemble.
- Stable 9 ns endpoints with repeat disagreement motivate genuinely different protein
  starting conformations.
- Stable, agreeing endpoints still cannot distinguish shared trapping from force-field or
  experimental-reference mismatch.
- No result from this pilot changes the frozen historical gate.

## First-light operational amendment — 2026-09-15, before production

The first real 10 ps continuation passed GROMACS's checkpoint/output checksum gate and
completed normally. GROMACS preserved every numeric XVG row through 3499 ps but regenerated
the checkpoint-boundary row at 3500 ps; the maximum reduced-potential difference in that
single column was 0.1774940004106611. The other 3500 historical rows were byte-for-byte
unchanged. This is the documented restart behavior in which output is restored to the last
checkpoint before continuation.

The integrity rule is therefore narrowed, not relaxed generally: source validation still
requires exact 500–3500 ps agreement before execution; after continuation, 500–3499 ps must
remain exact, only the 3500 ps boundary may differ, and its maximum absolute reduced-potential
difference is recorded in every output NPZ. The regenerated boundary is never analyzed: the
combined matrix keeps the checksum-frozen source's 3500 ps column and admits continuation
columns beginning at exactly 3501 ps. Any earlier change, gap, duplicate, or non-finite value
remains fatal.
