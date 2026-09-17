# F64A selection-policy sensitivity — fixed 2026-09-16 before execution

The user authorized this CPU-only package after inspecting the
[disjoint-block results](f64a_block_analysis_results.md). It is exploratory and
post-result, not a preregistered convergence test. No new GPU sampling is authorized.

Completed on SCC: see the [preserved report and results](f64a_selection_sensitivity_results.md).
The design below records the analysis as fixed before that run.

Question: does final-half movement persist when within-block cutoffs are fixed,
rather than re-estimated separately for every window and half?

## Fixed comparison

Use all 60 validated completed production NPZs, unchanged. Compare raw columns
`[6001:7501]` against `[7501:9001]` in each of the three repeats. Retained clocks
are 6.001–7.500 and 7.501–9.000 ns; no endpoint is duplicated.

Three policies are fixed in
[`cpu_sensitivity.yaml`](../pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml):

1. Historical adaptive trim-and-thin, unchanged.
2. Zero additional trimming within each half, with correlation thinning.
3. Discard the first 25% (375 columns) of each half, with correlation thinning.

For fixed policies, estimate statistical inefficiency separately for sampled-state
reduced potential and the same historical neighbor energy difference, **on the suffix
after the fixed cutoff**. Use the larger g and the historical subsampling function.
Do not run adaptive equilibration detection for either fixed policy. A failed or
non-finite correlation estimate aborts the whole report; do not silently replace it.
The 25% cutoff is a declared sensitivity choice, not proof of equilibration.

Eighteen MBAR fits record the same conditional uncertainties, aggregate/local
hysteresis, overlap matrices and window selection diagnostics as the previous report.
Also record forward/reverse adjacent Zwanzig estimates (both oriented lower→upper)
and their paired half-to-half changes. These estimates can be noisy; local hysteresis
or adjacent changes do not decompose the MBAR endpoint drift or identify its cause.

Report repeat means/SEM/ranges, propagated MBAR error and larger-of uncertainty,
and paired movement/SEM separately for each policy. Policies reuse the same data;
do not count them as independent replications or attach significance tests. All
uncertainties remain conditional on sampled modes, not model/reference bias.

## Integrity and execution

Require SCC primary and block reports to match their preserved repository references.
Require matching simulation identities, protocol and successful prior reproduction.
The six adaptive fits must reproduce all six scalar diagnostics within 1e-8 kcal/mol
and exact per-window selection records/counts. This numerical tolerance is not a
scientific convergence threshold. Hash the 60 NPZs, manifest, cluster reports and
frozen references before/after execution; also require every production NPZ hash to
match its recorded input hash in the frozen block report. Capture clean committed analyzer/config/
reference/shared-source identity and runtime versions separately from simulation identity.
Failure aborts without a partial result. Atomic publication refuses existing output.

After committing/transferring the package, and only with the GPU array finished:

```bash
cd /projectnb/rise-batteries/bode/rise-project
source scripts/scc_env.sh
python -m src.analysis.f64a_sensitivity \
  --config pilots/f64a_sampling_extension_v1/cpu_sensitivity.yaml
```

Output: `results/pilots/f64a_sampling_extension_v1/production/analysis/f64a_selection_sensitivity_v1.json`.
Do not rerun or overwrite either prior report. The analyzer deliberately refuses
execution from uncommitted relevant paths. Snakemake equivalent:

```bash
snakemake f64a_extension_sensitivity --snakefile workflow/Snakefile --cores 1 \
  --allowed-rules f64a_extension_sensitivity
```

The local rule consumes completed inputs only; missing matrices or block reports
are fatal and cannot schedule GPU producers. Neither GPU pilot nor historical
pipeline config is changed. Persistent downward movement across policies would
strengthen the case for a separately designed, authorized bounded continuation;
strong dependence on selection would motivate resolving estimator sensitivity first.
No policy is promoted to the correct estimate merely because it appears stable.

## Local package verification

Before transfer, `python -m pytest -q tests/` passed all 184 tests. Ruff and
`git diff --check` passed. The default Snakemake dry-run passed; the targeted
synthetic-completed-input dry-run contains only the CPU sensitivity producer and
fails when its completed block report is removed. Synthetic tests also verify
suffix-only correlation estimation, no adaptive detection for fixed policies,
adaptive reproduction failures, prior input hash mismatches, mid-run mutations,
refusal to overwrite, and preservation of both prior reports. These tests do not
constitute real F64A sensitivity results; those require the SCC matrices.
