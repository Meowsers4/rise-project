# F64A disjoint-block CPU diagnostic — fixed 2026-09-16 before execution

## Status and question

The user authorized this CPU-only follow-up after viewing the completed 3/6/9 ns pilot
report. It is **exploratory**: the cumulative results are already known. This analysis
has now been run on SCC; the supplied report and interpretation are preserved in
[the block results](f64a_block_analysis_results.md). The real extended matrices
remain on SCC and have not been solved locally.
The [prior results](f64a_extension_results.md) and their supplied JSON are preserved,
not overwritten. This design selects no new GPU experiment and introduces no convergence
pass/fail threshold.

Question: does agreement improve between successively later disjoint blocks, including
within the last 3 ns, and how much data does the interval-specific trimming/thinning
retain? Per-window energy-space diagnostics may locate problematic ladder regions;
they do not identify a protein conformational mechanism.

## Fixed input and intervals

Input is all 60 validated production NPZs for the completed F64A folded apo-2SH extension.
Use the existing complete protocol/provenance/code/config/shape validator. Read the
stage manifest and original SCC analysis; require the latter to agree with the preserved
report. Record SHA-256 checksums of the 60 NPZs, stage manifest, frozen report, new CPU
config, analyzer source, shared estimator source, and the analyzer's Git/runtime identity.
Keep simulation identity separate from later analysis identity. GPU pilot config and
base pipeline config remain unchanged.

Raw retained-column index i represents original simulation time 500+i ps. Slices are
half-open; no endpoint is duplicated across a disjoint comparison:

| Label | Raw-column slice | Retained clock, ps | Count/window |
|---|---|---|---:|
| 0–3 ns | `[0:3001]` | 0–3000 | 3001 |
| 3–6 ns | `[3001:6001]` | 3001–6000 | 3000 |
| 6–9 ns | `[6001:9001]` | 6001–9000 | 3000 |
| First half of final block | `[6001:7501]` | 6001–7500 | 1500 |
| Second half of final block | `[7501:9001]` | 7501–9000 | 1500 |

The final 3 ns block and its halves overlap by design and are **not** three independent
replications. Each configured contrast uses genuinely disjoint raw-column sets within
one repeat, but time separation does not guarantee statistical independence of slow modes.

## Estimation and records

Use the same neighbor-aware adaptive equilibration trimming and correlation thinning as
the primary report: estimate t0 and g from the sampled-state reduced potential and its
neighbor energy difference, take the larger t0 and g, then subsample. Instrument the
existing shared implementation without changing its selection algorithm.

For five intervals × three repeats (15 MBAR solves), report folded ΔG and conditional
MBAR uncertainty, signed-net/sum-absolute/maximum-local hysteresis, minimum adjacent
overlap, complete overlap matrix, and solver notes. Do not omit failed fits; abort without
writing a plausible partial result. Preserve 20 per-window diagnostic records for each
fit (300 total): both t0 values and g estimates, chosen cutoff and g, retained count,
first/last selected raw columns and corresponding times. Counts are heuristic retained
samples, not a proof of independent conformational observations. Record 19 signed local
adjacent discrepancies per fit (285 total).

Summarize repeat mean, SEM, range, propagated MBAR uncertainty, and their larger-of rule
for each block. Report signed/absolute paired movements for first→middle, middle→final,
and final-half→final-half (nine repeat-level contrasts). Summarize their paired mean and
SEM descriptively. Do not attach independent-sample p-values or promote three repeats
into a broad diagnostic-validation claim. If late estimates change while cutoffs also
change materially, record that estimator-selection sensitivity remains unresolved;
do not attribute all movement to a precisely timed conformational transition.

The 0–3 and 6–9 ns results must reproduce the corresponding primary-report estimates,
uncertainties, diagnostics, and selected sample counts. The config's 1e-8 kcal/mol
tolerance is numerical reproduction tolerance, **not** a scientific convergence threshold.

## Execution and output

Run only after the production array is finished; confirm `qstat` is empty for job 7589742
before updating the checkout. Use `source scripts/scc_env.sh`, then:

```bash
python -m src.analysis.f64a_blocks \
  --config pilots/f64a_sampling_extension_v1/cpu_blocks.yaml
```

This is a CPU-only command: no `qsub`, GROMACS, structural rebuilding, or new sampling.
It validates existing inputs and refuses to overwrite the output
`results/pilots/f64a_sampling_extension_v1/production/analysis/f64a_blocks_v1.json`.
Snakemake rule `f64a_extension_blocks` represents the CPU dependency graph without
making changed CPU code a reason to rerun any GPU job; missing completed inputs are fatal.
When using Snakemake, explicitly restrict permitted producers to this CPU rule:

```bash
snakemake f64a_extension_blocks --snakefile workflow/Snakefile --cores 1 \
  --allowed-rules f64a_extension_blocks
```

Interpret the latest disjoint comparison and per-window records before designing any
further GPU experiment. Persistent late movement may motivate a new bounded sampling
test; stable blocks with persistent repeat disagreement may motivate genuinely different
protein starts. Stable, agreeing blocks still cannot exclude shared trapping or resolve
model/reference mismatch. Preserve the historical failed gate and folded-only scope.
