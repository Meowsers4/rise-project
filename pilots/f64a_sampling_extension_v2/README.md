# Isolated F64A 9→15 ns continuation

User selected this bounded arm on 2026-09-17. Read
[PREREGISTRATION.md](PREREGISTRATION.md) before execution. This package never modifies
v1 or the historical folded/unfolded results. Mentor review is deferred.
See [VERIFICATION.md](VERIFICATION.md) for completed local checks and pending SCC work.

## 1. Deployment and read-only checks

Commit/push the package, then `git pull --ff-only` on SCC **only with no array in
flight**. All source matrices/checkpoints remain on SCC. From the repository root:

```bash
cd /projectnb/rise-batteries/bode/rise-project
source scripts/scc_env.sh
mkdir -p logs/f64a_continuation
qstat -u "$USER"
qacct -j 7589742 > logs/f64a_continuation/v1_qacct.txt
pquota > logs/f64a_continuation/quota.txt
acctool -project rise-batteries -detail -1 -balance 09/01/26 yesterday \
  > logs/f64a_continuation/balance.txt
python -m src.fep.f64a_continuation verify-tools
python -m src.fep.f64a_continuation validate-source
```

`validate-source` is CPU-only/read-only: 60 NPZs, 360 append/TPR/stamp files and
four ancestry files, with exact histories and actual TPR/checkpoint inspection.
Missing files, truncated energies, mismatched reports or unexpected boundaries stop
before staging. No new simulation or destructive cleanup is performed by these commands.

## 2. Resource admission — fill actual reviewed values, not placeholders

Review measured v1 accounting and project quota, confirm the current SU/core-hour
rate for the L40S placement, and select aggregate caps and finite attempts. No values
below have defaults; do not execute the example literally:

```text
python -m src.fep.f64a_continuation admit \
  --qacct-file logs/f64a_continuation/v1_qacct.txt \
  --quota-file logs/f64a_continuation/quota.txt \
  --balance-file logs/f64a_continuation/balance.txt \
  --available-gib REVIEWED_PROJECT_QUOTA_HEADROOM_GIB \
  --su-balance CURRENT_PROJECT_BALANCE \
  --su-rate CONFIRMED_SU_PER_REQUESTED_CORE_HOUR \
  --gpu-hour-cap APPROVED_TOTAL_GPU_HOURS \
  --su-cap APPROVED_TOTAL_SU \
  --max-attempts APPROVED_FINITE_ATTEMPTS_PER_TASK
```

This captures raw evidence and source hashes in immutable `admission.json`, estimates
storage requirements with the declared margin, and derives scheduler walltimes from
the aggregate caps. It fails if measured runtime/headroom, quota or SU budget cannot
fit. GiB means bytes/1024³; convert quota display units accurately. Input values are
explicitly human-reviewed, not automatically extracted from ambiguous quota tables.

The requested-core-hour reservation uses the confirmed node rate, not `qacct` CPU
utilization. See [BU project accounting](https://www.bu.edu/tech/support/research/account-management/manage-project/).

The historical 36 GPU-hour estimate is not permission to spend that amount or a
substitute for actual accounting. Do not edit code/config after admission. If a
preflight fails, preserve it and diagnose rather than overwriting/removing evidence.

## 3. First light only

```bash
python -m src.fep.f64a_continuation stage --first-light
python -m src.fep.f64a_continuation submit --first-light
qstat -u "$USER"
```

Do **not** call `qsub` directly: only the guarded submit command records the job ID
and overrides `h_rt` using the admission. Submission holds the job until its identity
is recorded, then releases it. If recording/release fails, inspect the reported held
job ID and logs; do not duplicate it. After this job finishes:

```bash
tail -n 100 logs/f64a_continuation/f64a_v2_light.*
python -m src.fep.f64a_continuation check-first-light
```

Expected: `(20, 9011)`, finite, exact frozen 9001-column source prefix. A technical
failure requires inspection; never launch production just because the job left qstat.
No pull/edit while first light or production is running. An existing submission
record refuses duplication; use the recorded job ID to inspect accounting first.

## 4. Production, after first light passes and resource snapshot is reviewed

```bash
python -m src.fep.f64a_continuation stage
python -m src.fep.f64a_continuation submit
qstat -u "$USER"
find results/pilots/f64a_sampling_extension_v2/production/fep/F64A/folded \
  -maxdepth 1 -type f -name 'w*_r*.npz' | wc -l
```

Wait for the array to finish, confirm 60 outputs and inspect qacct/log failures. A
finite per-task attempt ledger prevents repeated reschedules from exceeding admission.
Do not automatically submit another array or delete a checkpoint on failure.

After a technical failure has been diagnosed and all prior jobs in that mode have
finished, explicit remaining-budget retries use the same guarded submit command:
`python -m src.fep.f64a_continuation submit --retry-tasks 5` (an example task ID,
not an instruction to retry them). For first light use `submit --first-light
--retry-tasks 1`. It refuses active duplicates, existing outputs and exhausted
attempts, preserves prior submission records, and resumes the private checkpoint.
Multiple retry IDs must form a contiguous range (for example `5,6`), following
SCC's documented `-t start-end` syntax; do not include successful tasks in that range.

## 5. Completed-input CPU analysis

```bash
python -m src.fep.f64a_continuation validate-outputs
python -m src.analysis.f64a_continuation
```

Equivalent CPU-only workflow target:

```bash
snakemake f64a_v2_analysis --snakefile workflow/Snakefile --cores 1 \
  --allowed-rules f64a_v2_analysis
```

Missing completed inputs are fatal, never GPU producers. Output:
`results/pilots/f64a_sampling_extension_v2/production/analysis/f64a_continuation_v2.json`.
All 54 fits and v1 reproduction must succeed; existing reports cannot be overwritten.
The stage/GPU rules represent dependencies in Snakemake; actual GPU dispatch remains
the guarded SGE submission above, not local execution of the GPU target.

Stop at 15 ns. Compare all policies and repeats without choosing the result that
looks stable. This is folded-leg time dependence, not a folding ΔΔG or gate rerun.
