# Technical pre-execution amendment — 2026-09-17

Status: v2 source preflight stopped before resource admission, staging or submission.
The original committed text fingerprint excluded only `inputrec.nsteps`. It also
included a non-active lambda scalar in the TPR header, causing a false mismatch.
No energy history, sampling protocol, analysis policy or endpoint is amended.

## Observed evidence

After successful w0/r0–r2 checkpoint checks, the w1/r0 model comparison failed.
The user ran GROMACS 2025.3 `check -s1 .../source/.../prod.tpr -s2 .../fep/.../extension.tpr
-tol 0 -abstol 0`. Its only reported difference was `nsteps` (1750000→4750000).
Topology, force-field parameters, box, coordinates and velocities had no reported
differences. A dump diff with `nsteps` removed showed only:

```diff
header:
-   lambda = 5.880000e-02
+   lambda = 0.000000e+00
```

This explains why window zero passed: its original header lambda was already zero.
These are user-supplied diagnostics of existing TPRs, not new simulation results.

## Source-based explanation

GROMACS 2025.3 [`convert_tpr.cpp`](https://github.com/gromacs/gromacs/blob/v2025.3/src/gromacs/tools/convert_tpr.cpp)
creates a state, reads the TPR, increments the input record's step count, and writes
the state back. [`tpxio.cpp`](https://github.com/gromacs/gromacs/blob/v2025.3/src/gromacs/fileio/tpxio.cpp)
reads the header lambda separately, without copying it into the state's lambda
array; writing populates that header scalar from the state. The constructor in
[`state.cpp`](https://github.com/gromacs/gromacs/blob/v2025.3/src/gromacs/mdtypes/state.cpp)
initializes that array to zero. `initialize_lambdas` in the same file sets active
lambda components from the FEP input record. Thus the observed header reset is
consistent with the conversion path, not a change to the input-record schedule.
The pilot still resumes the original checkpoint with GROMACS append checksums.

## Narrow correction and retained checks

- Exclude `nsteps` only inside the `inputrec` section.
- Normalize only the exact scalar named `lambda` inside the top-level `header`
  section. Never normalize lambda/state/schedule fields in the input record,
  topology, coordinates, velocities or checkpoint.
- Require unique input-record/header/topology sections and one finite header
  lambda in [0,1]. Converted TPRs must carry the observed zero header value.
- Record original and converted header values, the explicit normalization policy,
  and a clearly named `model_sha256` in the admission's TPR inspection records.
  All other dump fields must still hash identically.
- Keep the active state-index check, 1 ps energy cadence, final TPR/checkpoint
  times, no-coordinate retention, full source hashes and exact energy histories.

Regression tests reproduce the supplied header difference and prove that changes
to active lambda settings/schedule, topology, atom count, coordinate/velocity data
and other header fields still change the fingerprint. Invalid header metadata fails.

Commit/deploy the correction only with no array in flight, then rerun all 60 source
checks. This amendment does not waive the remaining preflight or first-light gates.
Do not remove or overwrite any existing approval/stage evidence if one exists.
