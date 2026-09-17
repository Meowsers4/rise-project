# Preparation and verification — 2026-09-17

The user selected the isolated 9→15 ns arm. Package, registration, guarded submission,
finite resource/attempt contract and completed-input CPU workflow are implemented.
No new GPU job or scientific result has been produced by this preparation.

Local checks in the working repository:

- `python -m pytest -q tests/`: **228 passed**, including synthetic continuation,
  restart-boundary, source-preservation, budget, submission and 54-fit report tests.
- Ruff on both new modules and their tests (plus the two updated workflow fixtures): passed.
- `bash -n` on both new SGE scripts: passed.
- `git diff --check`: passed.
- Default Snakemake dry run: passed. The v2 CPU target also dry-runs against
  synthetic completed inputs without any GPU rule; missing inputs fail instead
  of scheduling simulation. These are workflow tests, not real energy estimates.
- Existing v1 code/config hashes match the preserved sensitivity report.

SCC checks still pending:

- Current job/accounting, project quota, balance and applicable node charging rate.
- All 60 actual source histories, archived/continuation TPR models and checkpoints.
- Explicit aggregate GPU-hour/SU caps and finite attempt allowance, then admission.
- Real 10 ps first light, before any production staging/submission.
- Production completion and the real fixed 54-fit CPU analysis.

An SSH BatchMode check could not authenticate to SCC from this session. No credentials
were requested or changed. No source staging, cluster resource admission or submission
was performed. Follow the [README](README.md) from the SCC shell after deployment;
review the evidence before choosing budget values. Mentor review remains deferred.
