# Pre-registration — G93A disulfide diagnostic (SS vs 2SH)

**Status: PRE-REGISTERED. Signed off by the user 2026-09-12; config staged on branch
`diag/g93a-ss` (never to be merged). No SS window has been run.** `CLAUDE.md` rule 1 required
that sign-off because v1 is defined as the disulfide-reduced form. Everything below — baseline,
endpoint, and all five outcome readings — was written and committed (`6d3926f`, `5430dac`)
before any SS window existed, so the interpretation cannot be adjusted to the result.

> **The automated checks cannot see this experiment.** The full test suite passes unchanged
> (105 passed, 12 skipped) with `keep_disulfide_reduced: false`, and the protocol fingerprint is
> byte-identical to the 2SH baseline (§5). Nothing in the repo will tell you which physical
> state a window came from. The branch discipline (§7 step 2) and the MANIFEST (§7) are the
> only provenance this arm has.

## 1. The question

[`reference_state_audit.md`](reference_state_audit.md) established that every gate control was
measured with the Cys57–Cys146 disulfide **intact**, while every calculation was run on the
**reduced** form. Whether that mismatch costs anything in a ΔΔG — a difference of differences,
where a WT/mutant-identical disulfide contribution cancels — is currently an assumption in both
directions. This measures it once, on one variant.

**It does not attempt to fix the gate.** Matching the disulfide could repair at most 5 of the 7
usable gate points: F64A is excluded by a pre-registered convergence criterion, and G93S/G93V
are DSC-on-apo-dimer, a different observable that no redox change addresses. The gate failed
and pivoted (r = 0.326 < 0.60); README §10 commits the project to a methods/limitations result.

## 2. Design — one variable

| | 2SH baseline (exists) | SS diagnostic (proposed) |
|---|---|---|
| variant | G93A | G93A |
| protocol hash | `822108e9db71124d` | **`822108e9db71124d`** (verified identical, §5) |
| folded sampling | 3 ns + 0.5 ns equil | same |
| λ-windows / replicates | 20 / 3, independent boxes | same |
| **Cys57–Cys146** | **reduced** | **intact** |

G93A because it has the cleanest cycle closure in the project (0.045), a Lindberg reference
verified at source (ΔG 3.03 − 0.60 = 2.43, matching `variants.csv`), and because `CLAUDE.md`
already uses it as the "converged but wrong by 1.17 kcal/mol" exemplar. If the disulfide moves
that error, this is where it should be visible.

## 3. Frozen baseline

Re-derived from archived windows 2026-09-12 ([`stop_rule_reanalysis.md`](stop_rule_reanalysis.md)),
reproducible in three environments. **These numbers are fixed as of this document.**

| | value |
|---|---|
| ΔΔG | 1.2587 ± 0.0612 (exp 2.43, **error −1.171**) |
| per-replicate ΔΔG | 1.3782 / 1.1759 / 1.2219 (spread 0.2023) |
| **folded ΔG (primary endpoint)** | r0 **+9.9138**, r1 **+9.6637**, r2 **+9.7211** — **mean +9.7662, sd 0.1311** |
| unfolded ΔG | +8.5356 / +8.4877 / +8.4992 |
| cycle closure | 0.0450 |
| min adjacent overlap | 0.02486 |

**Primary endpoint is the folded-leg ΔG, not ΔΔG.** The unfolded reference is a capped
tripeptide (ACE-G92-G93-A94-NME) containing no cysteine, so it is disulfide-free in both arms by
construction and cancels exactly. Comparing folded ΔG isolates the disulfide; comparing ΔΔG
re-introduces unfolded-leg sampling noise (~0.1 kcal/mol) for nothing. ΔΔG is reported as a
secondary, human-readable figure.

## 4. Pre-declared interpretation

Let **Δ = mean folded ΔG(SS) − 9.7662**, compared against the 2SH replicate sd of 0.1311.

| outcome | reading | consequence |
|---|---|---|
| **\|Δ\| < 0.3** | The disulfide does not materially change this mutation's folded-state cost. | The reference-state mismatch is **documented but not quantitatively important for ΔΔG at this site**. LiveCoMS keeps the mismatch as a finding about benchmark construction, and adds a direct negative test. Second candidate mechanism excluded — same value as the sampling test. |
| **0.3 ≤ \|Δ\| < 1.0** | Real but partial. | Report as a bounded systematic; state that it does not close the 1.17 gap alone. No re-gate. |
| **\|Δ\| ≥ 1.0, sign reduces the error** | The mismatch is a major contributor. | Strongest outcome. Still **no re-gate** — one variant cannot re-open a pre-registered gate. Becomes the argument for a properly pre-registered state-matched campaign as future work. |
| **\|Δ\| ≥ 1.0, sign increases the error** | Real and adverse. | Report honestly. Strengthens the limitations note; weakens any future case for state-matching. |
| **SS run fails to converge** (closure > 1.0, or overlap collapse) | The SS folded state is a harder sampling problem. | Report as a null with the diagnostics; do not retry at longer sampling without a new pre-registration. |

**No outcome licenses re-running the gate, re-weighting a control, or removing one.**

## 5. Verified before writing this

Computed locally 2026-09-12 against the committed engine:

- With `ns_per_window.folded: 3` and `equilibration_ns.folded: 0.5`, the folded fingerprint is **`822108e9db71124d`** — byte-identical to the G93A archive. Reverting those two values is necessary **and sufficient** for protocol comparability.
- Flipping `keep_disulfide_reduced` **does not change the fingerprint**: the disulfide is a topology property, absent from the `.mdp`, and `protocol_extra()` carries only `independent_replicate_systems`. **The hash cannot witness this experiment** — §6 trap 2.
- `_pdb2gmx_stdin()` ([pmx_engine.py:294](../src/fep/pmx_engine.py#L294)) returns `""` when the flag is false, so pdb2gmx falls back to its default and forms C57–C146 by SG–SG distance. The disulfide-free guard at [line 574](../src/fep/pmx_engine.py#L574) is gated by the same flag.
- No gate variant's tripeptide contains Cys57 or Cys146 (closest: I149A at 148–150), confirming §3's cancellation argument.

## 6. Three traps that would make this measure nothing

**1. The cached system silently defeats the flag.** `build_system` fast-paths on the
`SYSTEM_READY` marker ([pmx_engine.py:508](../src/fep/pmx_engine.py#L508)), which encodes no
topology or config. `results/fep/G93A/folded/system_r{0,1,2}/` already exist **in the reduced
state**; flipping the flag and resubmitting would reuse those boxes and never re-run pdb2gmx.
*Mitigation: move the entire `results/fep/G93A` tree aside before submitting (§7).*

**2. The resume guard cannot protect this run.** `assert_resumable` compares protocol hashes —
and SS and 2SH **share** `822108e9db71124d`. A stale 2SH `prod.cpt` would therefore pass the
guard and be resumed under the SS config. This is precisely the 2026-08-11 G93A failure the
guard's own docstring describes, and the guard is blind to it here. *Mitigation: same as trap 1
— no stale run dir may exist. This is why the tree is moved, not merged.*

**4. `structure.disulfide` must stay `reduced` — do not "fix" it.** This is the most
counter-intuitive part of the design and the one most likely to be broken by a well-meaning
edit. `prepare_variant` **hard-raises** unless `structure.form == "apo"` and
`structure.disulfide == "reduced"` ([build.py:187](../src/prep/build.py#L187)), so setting it to
`oxidized` does not produce the SS state — it aborts the build.

The SS state is produced by a different route entirely. Stage 1 calls `strip_disulfide_bonds`,
which removes the SG–SG **bond from the OpenMM topology only** and does not touch coordinates
([build.py:117](../src/prep/build.py#L117)) — the two SG atoms stay at their crystal separation
of ~2 Å. pdb2gmx then re-detects disulfides **by SG–SG distance**, and with the `-ss` answers
suppressed (`keep_disulfide_reduced: false` → `_pdb2gmx_stdin` returns `""`) it re-forms
C57–C146 on its own. The SS arm is therefore the *inverse of the guard*: it works by declining
to answer a prompt, not by declaring an oxidized state anywhere in config.

Consequence: `structure.disulfide: reduced` remains correct and unchanged throughout this
experiment, and the only config key touching redox is `fep.keep_disulfide_reduced`.

**3. Sampling must be reverted too, or it confounds.** The committed config is 9 ns / 2.0 ns
folded. Running at 9 ns would vary disulfide **and** sampling simultaneously and break
comparability with the 1.17 baseline. Three config values change, not one.

## 7. Procedure

Three config values, one tree move, one smoke test, one array.

```bash
# ---- 0. on the SCC, from the repo root ------------------------------------------
cd /projectnb/rise-batteries/bode/rise-project
qstat -u bodeb                      # MUST be empty; never pull or reconfigure mid-array
git pull                            # picks up this document

# ---- 1. protect the 2SH baseline (traps 1 and 2) --------------------------------
# results/ is gitignored and this tree is the 2SH baseline. It is archived in two
# places already, but move rather than delete.
mv results/fep/G93A results/fep/G93A_2SH_baseline
find results/fep/G93A_2SH_baseline -name 'w*_r*.npz' | wc -l    # expect 120
ls results/fep/G93A 2>/dev/null && echo "STOP: G93A still present" || echo "clean"

# ---- 2. the three config values, ON A BRANCH THAT IS NEVER MERGED ---------------
#   fep.keep_disulfide_reduced : true -> false
#   fep.ns_per_window.folded   : 9    -> 3
#   fep.equilibration_ns.folded: 2.0  -> 0.5
#   fep.structure.disulfide    : UNCHANGED ("reduced") -- see trap 4
#
# CRITICAL: these must NOT land on main. If they do, every later array -- including the
# week-2/3 F64A and G93V runs -- silently executes SS at 3 ns. Use a throwaway branch:
#   (locally)  git checkout -b diag/g93a-ss && <edit> && git commit && git push -u origin diag/g93a-ss
#   (here)     git fetch && git checkout diag/g93a-ss
# Never merge it. After the run: git checkout main  (step 7).
# Do NOT hand-edit config on the SCC (CLAUDE.md).
git branch --show-current      # MUST print diag/g93a-ss before submitting

# ---- 3. smoke test FIRST: the SS path has never run end-to-end ------------------
# The system build is what is being tested (does pdb2gmx form the bond?), and that runs
# before any mdrun, so --smoke is enough and takes minutes rather than 15.
source scripts/scc_env.sh
python -m src.fep.window --variant G93A --leg folded --window 0 --rep 0 --smoke \
  --config config/pipeline.yaml --out "$TMPDIR/ss_smoke_w0_r0.npz"
#   npz goes to scratch, NOT into results/: run_window never skips an existing output, so a
#   smoke-hash window left in results/fep/G93A/folded/ would make _check_single_protocol
#   reject the whole variant later.

# The guard is DISABLED in this arm (keep_disulfide_reduced: false), so this grep is the
# ONLY verification that the topology is what we intend. Check the file the guard would
# have read, and check BOTH directions:
T=results/fep/G93A/folded/system_r0/hybrid.top
grep -n "CYS2\|CYX" "$T"          # expect the C57/C146 pair present -> the bond formed
grep -c "CYS2\|CYX" "$T"          # necessary but NOT sufficient on its own
#   Then confirm no SPURIOUS bond: SOD1 has four cysteines (6, 57, 111, 146). Only 57-146
#   may be bridged; C6 and C111 must remain free thiols with HG. If C6 or C111 appear as
#   CYS2/CYX, pdb2gmx has over-bonded and the run is invalid -- STOP.

# Clean up: --smoke stamps the run dir with a DIFFERENT protocol hash, so the array's real
# w0_r0 task would hit assert_resumable and refuse. Remove the run dir; keep the system dir
# (identical either way, and rebuilding costs time).
rm -rf results/fep/G93A/folded/w0_r0

# ---- 4. the array ---------------------------------------------------------------
mkdir -p logs/fep
qsub -v VARIANT=G93A scripts/submit_array.sh        # 120 tasks, ~18 GPU-h

# ---- 5. monitor ------------------------------------------------------------------
echo "$(date +%H:%M) done: $(find results/fep/G93A -name 'w*_r*.npz' | wc -l)/120 | live: $(find results/fep/G93A -name prod.log -newermt '-2 minutes' | wc -l) | queue: $(qstat -u bodeb | grep -c sod1_fep)"

# ---- 6. analyse into a SEPARATE file; never overwrite the baseline --------------
python -m src.fep.analyze --variant G93A --config config/pipeline.yaml \
  --out results/fep/G93A/ddg_SS.json

# ---- 7. label the tree, restore the baseline, leave the branch ------------------
mv results/fep/G93A results/fep/G93A_SS_diagnostic
mv results/fep/G93A_2SH_baseline results/fep/G93A
git checkout main                   # diag/g93a-ss must never be merged
```

**The manifest matters more than usual here**, because the hash cannot witness the state
change (§5). `results/fep/G93A_SS_diagnostic/MANIFEST.md` must record: disulfide **intact**,
`keep_disulfide_reduced: false`, the commit SHA, the date, and that its hash
`822108e9db71124d` is **shared with the 2SH baseline and therefore proves nothing about the
redox state**. Without that note, a future `analyze` run could merge the two sets and the guard
would not object.

It must also record that **both legs carry `822108e9db71124d`** — the folded leg is physically
SS, the unfolded one is state-free — so the shared hash across legs is expected here and is not
evidence of contamination. And it must state that `structure.disulfide` was `reduced`
throughout (trap 4), so the manifest does not look self-contradictory to someone who has not
read this document.

## 8. Cost

| | |
|---|---|
| folded | 60 windows × ~15 min = **~15 GPU-h** |
| unfolded | 60 × ~2–3 min = **~3 GPU-h** |
| **total** | **~18 GPU-h**, one array, well inside `h_rt=12:00:00` per task |

**Both legs run; there is no folded-only mode.** `submit_array.sh` hardcodes `#$ -t 1-120` and
decodes every task across `fep.legs`, and the runtime guard recomputes `legs*windows*reps = 120`
and refuses `-t 1-60`. **Do not hand-edit `-t` without also fixing that guard** — they are
coupled deliberately, and a mismatch exits 2 on every task. So the array is 120 windows, not 60.

The unfolded leg is re-run rather than reused. It is physically identical (no cysteine in the
tripeptide), so reuse would be defensible — but it costs 3 GPU-h to avoid copying result files
between trees, and copying is how provenance mistakes start. The primary endpoint ignores it
either way.

## 9. What this is not

- Not a gate point, not a gate re-run, not an input to `evaluate_gate`. `validation.gate_subset` is unchanged.
- Not a switch of the project's state. v1 remains apo-2SH; this is a single diagnostic arm.
- Not a holo/metal-bound run — the thing `CLAUDE.md` rule 1 actually forbids.
- Not a claim about the other seven variants. One site, one mutation, one force field.
