# Direction review after the gate pivot — 2026-09-11

Summary of an external review of the project's post-pivot options, received 2026-09-11
against `main` at `85c69bd`. It answers two questions: (1) is the planned paper supported
by the present evidence, and (2) if not, what replaces it.

**Headline: do not spend ~350 GPU-hours on gate attempt 2.** There is a defensible
one-month archival submission, but only as a bounded methods-and-limitations note. The
present evidence supports neither a SOD1 biology paper nor a general result about
convergence diagnostics.

---

## 1. What was and was not verified

| | |
|---|---|
| Reviewer ran | full test suite (117 cases from 106 functions across 8 files, 37 in `test_fep.py`) and `snakemake -n` (978 planned jobs, 960 FEP windows). Both passed. No arrays submitted; tracked repo unchanged. |
| Reviewer did **not** | access the SCC. Local `results/` holds only `figures/`. |
| Therefore | **every reported SOD1 result in this repo remains unverified against raw output.** The review reproduced arithmetic from committed documentation, not MBAR estimates from NPZ/XVG files. |
| Re-checked locally 2026-09-11 | config values, the stale `0.407` comment, `analyze.py` uncertainty, fingerprint scope, `src/md/run.py`, Snakefile rule bodies, all gate arithmetic below. All confirmed. |
| Not re-derived here | protocol hashes and the test/DAG counts — no local env has `yaml`/`numpy`. Those figures are the reviewer's. |

## 2. Corrections to committed documentation

| Documented claim | What the evidence establishes |
|---|---|
| Gate failed at 0.407 **or** 0.326 | The seven rounded rows give r = 0.32597, RMSE = 2.12355, SSE = 31.5663. The `0.407` in [`pipeline.yaml:167`](../config/pipeline.yaml) is **stale**. |
| Protocol is 18 windows / 5 replicates, hash `822108e9db71124d` throughout | Config specifies **20 windows / 3 replicates**. Per-leg sampling split the hashes: folded `cf1e632168579261`, unfolded `822108e9db71124d`. |
| Snakemake rule bodies are all TODO | **False.** TODO comments remain, but rules invoke real modules and the gate DAG constructs. Genuinely missing: prescreen execution and `src/md/run.py`. |
| F64A's 1.02 replicate range vs ±0.34 proves understated uncertainty | **False.** For 6.61 / 7.62 / 6.60, SEM = 0.338. [`analyze.py:374`](../src/fep/analyze.py#L374) reports `max(propagated MBAR error, replicate SEM)`. Range and uncertainty-of-the-mean are different quantities. |
| Provenance + hash machinery makes fabricated or mixed-protocol results impossible | **Overstated.** Engine provenance is a string check; [`protocol_fingerprint`](../src/fep/pmx_engine.py#L740) covers the generated MDP plus selected extras, not full topology/force-field provenance. Changing `pmx_forcefield` to `oplsaamut` left the fingerprint unchanged. |
| Wells's r = 0.81 is the comparable folding benchmark | Wrong observable and state: that number is apo-SS **dimerization** (0.81 AMBER alone, 0.92 AMBER/OPLSAA averaged). Its 250 ns limitation is about residue-level dynamical convergence. |
| F64A's experimental −0.20 is suspect on burial/atom-count grounds | The value is in Nordlund & Oliveberg's original Table 1, measured at pH 6.3 on a C6A/C111A/F50E/G51E background. That is a **reference-compatibility** issue to audit — not grounds to call the experiment wrong. The equivalent audit for G93S/G93V is not done. |
| README's "Baby 2025" citation | Resolves to Gupta, Sun & Levy, published April 2026; benchmarks nuclease and T4 lysozyme. Substantial existing competition. |

External-source corrections (Wells, Nordlund & Oliveberg, Gupta) are the reviewer's, from
full texts rather than search snippets; they are not independently re-checked here.

## 3. Question 1 — the paper

**Yes to a limitations note. No to the stronger paper as currently described.**

Target: LiveCoMS *Lessons Learned*, which explicitly accepts failed studies and
reproducibility investigations but requires lessons useful beyond the one project
(multiple systems desirable). It requires a **presubmission inquiry — send it in week 1.**
This is a defensible specialist target, not an acceptance prediction.

The central claim, at the strength the evidence supports:

> In the documented apo-2SH SOD1 equilibrium-window pmx/AMBER calculations, small
> within-ladder forward/reverse hysteresis coexisted with disagreement between
> independently solvated replicates and failure against the prespecified experimental
> gate, so that diagnostic alone did not justify prospective variant triage.

Until SCC reanalysis succeeds this is supported by committed documentation, not by an
independently reproduced result.

**Do not strengthen it** to "the standard diagnostic contains no information about
accuracy" — seven points cannot establish that. F64A does not identify which replicate
sampled the correct basin; two agreeing replicates are not ground truth. The
implementation is absolute forward/reverse Zwanzig hysteresis summed along one ladder,
then maximized over legs and replicates ([`analyze.py:263`](../src/fep/analyze.py#L263))
— **not** a closed thermodynamic network. That distinction must appear in the manuscript.

Withdraw these causal claims:

- Tight unfolded-leg replication shows precision under the tripeptide model, not adequacy of that reference.
- Independent solvent boxes are not independent sampling of protein conformational basins.
- The G93 ladder localizes a discrepancy to one site; it does not identify backbone entropy as the cause.
- The bundled protocol changes do not isolate the effect of endpoint spacing, independent boxes, caps, or electrostatics.

The convergence-vs-correctness principle itself is prior art, now including a 2026
multi-system paper (Kang et al.). The contribution here is the reproducible failure
record and operational lessons — not the principle.

### Predictable reviewer objections

| Objection | What a month closes | What remains open |
|---|---|---|
| "These are documentation tables, possibly reflecting implementation or reference mismatches." | Regenerate estimates from NPZ/XVG; inspect topology, redox state, numbering, units, caps, provenance; audit experimental constructs. | Missing archives cannot be reconstructed from prose; reference-state discrepancies cannot be repaired by deleting controls. |
| "Known limitation, too little data." | Executable diagnostic example, full failed-run accounting, bounded sampling-sensitivity experiment. | Seven mutations at five sites; one target/force field/estimator. |
| "You claim a sampling mechanism without observing it." | Quantify time dependence, replicate disagreement, estimator sensitivity from retained energy data. | No retrospective structure without trajectories; no separation of sampling vs force-field vs reference bias; no Wells replication. |

For a limitations note, absent charge-changing coverage is a stated scope boundary. For a
paper claiming C1–C4 it is a missing central result.

## 4. Why attempt 2 is not worth buying

From the recorded folded-window timing: 60 × (892.8/3600) × (11/3.5) = **46.77 GPU-h**,
plus 60 unfolded tasks at 2–3 min → **48.8–49.8 GPU-h per variant**, or **341–348 for
seven**. At two continuously occupied GPUs that is 7.1–7.3 days; at 50% utilization,
about two weeks.

Two unreconciled numbers, to settle from scheduler accounting before trusting any budget:

- The "~23 GPU-h per old variant" figure in `HANDOFF.md` does not reconcile with the window timings. If 23 is real occupied time, the extrapolation is closer to **67 h/variant**.
- `cluster.gpus_available: 8` contradicts the practical throughput of 2.

Gate arithmetic, corrected: current SSE 31.5663, of which G93V contributes 17.3889 and the
other six 14.1774. Holding those six fixed, G93V needs absolute error **≤ 1.254**
kcal/mol — not ~1.5 as recorded in [`gate_attempt_1.md`](gate_attempt_1.md). That is
conditional arithmetic, not a prediction. **The reason to decline attempt 2 is opportunity
cost and unresolved model/reference questions, not a computed success probability.**

## 5. The bounded four-week plan

Filenames below are deliverables, not existing files.

| Week | Work | Output | Compute |
|---|---|---|---|
| 1 | Retrieve all 8 gate variants and superseded archives; reproduce historical analyses under their matching configs; audit primary reference states; recalculate F64A per-leg estimates, local hysteresis and replicate SEM; send LiveCoMS inquiry. | ✅ [`raw_result_reconciliation.md`](raw_result_reconciliation.md), ✅ [`stop_rule_reanalysis.md`](stop_rule_reanalysis.md) (stop rule **passes**), ✅ [`reference_state_audit.md`](reference_state_audit.md) (first pass; 7 numbers still need library PDFs), ⬜ LiveCoMS inquiry. | **0 GPU-h** (CPU only) — spent 0 |
| 2 | *Conditional on week 1 reconciling:* run F64A under the committed 9 ns folded / 2 ns equilibration protocol; keep the existing 3/0.5 ns unfolded protocol, 20 windows, 3 independent boxes. | `results/fep/F64A/ddg.json`, `results/convergence/F64A.json`, `docs/f64a_sampling_sensitivity.md` | 49–50 nominal (possibly ~67); ~25–34 occupied h at 2 GPUs |
| 3 | Run G93V under the same protocol. Analyze both regardless of direction; keep historical and follow-up estimates separate. | matching FEP/convergence JSON, `docs/g93v_sampling_sensitivity.md` | same as F64A |
| 4 | Protocol-history table, diagnostic comparisons, limitations manuscript; make analysis reproducible from archived inputs; mentor reproduces a main table. | manuscript + executable supplement + archived data | 0 |

Submission interface already exists:

```bash
qsub -t 1-120 -tc 2 -v VARIANT=F64A scripts/submit_array.sh   # week 2
qsub -t 1-120 -tc 2 -v VARIANT=G93V scripts/submit_array.sh   # week 3
```

`submit_array.sh` already pins L40S, one GPU, `h_rt=12:00:00`. The extrapolated folded
task is ~47 min, well under the cap. Keep `fep.checkpoint_interval_min: 5`.

Before either submission: archive the old result and checkpoint directories intact and
**never resume their checkpoints under the new folded protocol**. Pull only between
arrays. Keep all frozen thresholds and the gate subset unchanged. Preserve
`cluster.trajectory_retention: estimates_only` — this plan makes no structural claim.

Total ≈ **100 nominal GPU-h** (sensitivity toward 134). No VUS, no charge corrections, no
enhanced sampling, no extra force fields, no Stage 4. **Do not combine the two new
variants with five old-protocol variants into a new gate.**

**Week-1 stop rule:** if provenance-consistent reanalysis does not reproduce the
documented F64A low-hysteresis / replicate-disagreement example, withdraw this
recommendation. Preserve any useful software failure report; do not spend GPU time
rebuilding a preferred narrative.

## 6. Question 2 — the stronger replacement

Alternatives for the student's effort, not four parallel workstreams.

| Rank | Candidate | Venue / timeline | Decision |
|---|---|---|---|
| 1 | Test whether diagnostics justify stopping/reporting decisions on unseen systems; publish a reusable audit harness | JCIM, 8–12 weeks (JCTC only with a real inferential advance) | **Choose this** |
| 2 | DMS/predictor failure analysis separating abundance, activity, folding stability | Protein Science, 6–10 weeks | Cheapest backup, largely anticipated |
| 3 | C1: validated charge-changing SOD1 calculations | JCIM/Protein Science, 3–6 months | Good science, wrong month |
| 4 | Enhanced sampling of the G93 discrepancy | JCTC, 3–6+ months | Least suitable now |

### Candidate 1 — decision-focused audit of public FEP campaigns

The question: *at fixed budget, which diagnostics predict disagreement with an independent
repeat, and do they also identify inaccurate experimental predictions on unseen systems?*
Those are two distinct endpoints — repeat agreement tests reproducibility under a model;
experimental agreement additionally tests the model and reference comparability.

Data availability is unusually good. The public OpenFE tables (~7.3 MB total) hold 1,145
transformation rows over 49 targets (PyMBAR4), 1,199 rows with experimental comparison and
failure annotations (PyMBAR3), and 7,193 leg/repeat rows of cumulative estimates — three
repeats, both legs, per-leg estimates, statistical errors, minimum overlap. Differing row
populations need explicit reconciliation. Raw archives add reduced potentials, sample
counts, replica-state indices, timing. **Start at 0 GPU-h**; download raw archives
selectively to `/projectnb` (the JACS PTP1B archive alone is 2.6 GB — SCC home is 10 GB).

The potentially new content is a validated decision procedure, not a plot collection:
separate single-estimate uncertainty, replicate-mean uncertainty and experimental
discrepancy; evaluate on held-out targets keeping chemical series together; use only
information available at the proposed stopping time (no full-trajectory decorrelation
leaking into early-stop decisions); report false reassurance and retained fraction per
threshold; beat baselines of simulation length, MBAR uncertainty and repeat SEM; keep
excluded and incomplete calculations in the accounting; make every decision traceable to
source file, analysis version and protocol metadata.

Not an uncontested niche: SAMPL benchmarks sampling efficiency, adaptive-allocation methods
exist, and a public 2026 repository already studies calibrated cycle-closure detection on
OpenFE repeats. The differentiator must be prospectively usable stopping/reporting
decisions and their transfer across targets. What makes it ours is the operational failure
experience — wrong residue mapping, stale checkpoints, protocol mixing, optimistic
uncertainty, unsupported reference matching — but the harness needs complete input/topology
manifests, analysis provenance and immutable archived outputs. The current hash is useful,
not sufficient.

**Week-1 killer test:** on one development and one held-out system, reproduce published
terminal estimates and determine whether any available diagnostic improves prediction of
independent-repeat discrepancy beyond reported uncertainty and sampling length. If the data
cannot support that comparison, or the benefit vanishes held-out, stop before building
infrastructure.

Ruled out this month: regenerating barnase-109, nuclease or T4 lysozyme from scratch.
Published final predictions are insufficient for a hysteresis/overlap benchmark — per-run
diagnostics or raw work data are required. At SOD1 timing scale even 62 mutations imply
~3,100 GPU-h (an illustrative scale, not a timing prediction).

### Candidate 2 — notes

Joining the public Axakova scores using `project.mature_offset`: abundance scores exist for
all 38 panel VUS; activity for 37/38 (I151S absent); each file covers 91/92 across the full
panel. That is score availability, not confident classification — and Axakova already
compares the maps against DDGun and ML predictors including regional discordance. A new
correlation matrix is not novelty. Neither DMS abundance nor activity is experimental
folding ΔΔG, and **a failed FEP cannot adjudicate their disagreements.** Week 1 should
reproduce Axakova's existing analysis and identify what is left unanswered; absent a clear
residual question, stop. C4 is therefore ruled in as data work, ruled out as a standalone
FEP-concordance claim while the gate fails.

### Candidates 3 and 4 — why not now

**C1:** the 17 charge-changing controls inherit the same experimental-reference audit
problem, and "no prior charged-SOD1 FEP" is absence of evidence, not exclusivity. At
current sampling, 17 × ~50 ≈ **850 GPU-h** before box-size checks and correction
comparisons. Week 1 would establish analytical/topological correctness on a small charged
model and the feasibility of a separately authorized control sub-gate. **Do not bypass the
submission guard** on non-gate variants; a changed arm needs an explicit scoped design.

**C4/enhanced sampling:** the configured SCC GROMACS build has no MPI, and replica exchange
requires communicating simulations — the independent-window array architecture cannot
acquire HREX via an MDP flag. Four tempering states at the current folded budget would cost
~187 folded GPU-h per variant before overhead. The deeper risk is that backbone entropy may
not explain the discrepancy and an inadequate tripeptide reference remains possible.

**Do not** train an ML surrogate on the current FEP results — seven correlated,
inadequately validated labels would reproduce their deficiencies. Release the
provenance/pre-registration harness as an artifact within candidate 1 or the LiveCoMS note;
its existence and test count are not themselves a scientific advance.

**Reviewer's allocation:** one week to reconcile and preserve the SOD1 evidence, then
commit to candidate 1. If a submission this term is non-negotiable, follow the bounded
LiveCoMS plan instead. Attempting both alongside school recreates the same resource problem
in another form.

---

## 7. Conflicts with committed rules — for the user to resolve

Per `CLAUDE.md` ("surface the conflict instead of silently choosing"), **nothing in
`CLAUDE.md` or `HANDOFF.md` was edited.** The review contradicts these standing statements:

1. `CLAUDE.md` — "`ddg_err` has understated the true uncertainty 3–10x on every variant measured so far." The F64A instance of that inference compares a range to an SEM; `analyze.py` already takes the max of MBAR error and replicate SEM.
2. `HANDOFF.md` §7 — "the machinery is sound and the folded leg is not converged" and the F64A claim as "publishable whether or not the gate passes." The review supports the narrower statement in §3 above, and only after raw-output reanalysis.
3. `HANDOFF.md` §1 — the single protocol hash `822108e9db71124d` is now per-leg.
4. `HANDOFF.md` §1 — "~23 GPU-hours" per variant does not reconcile with the recorded window timings.
5. `HANDOFF.md` §5.1 / `gate_attempt_1.md` — the F64A/G93V/G93S "suspect experimental value" framing should become a reference-compatibility audit.
6. `gate_attempt_1.md` — "G93V must fall to about 1.5" should read **≤ 1.254**.
7. `pipeline.yaml:167` — the `r=0.407` comment is stale; the recorded failure is 0.326.
8. `README.md` — the "Baby 2025" citation resolves to Gupta, Sun & Levy (April 2026).

Items 3, 6 and 7 are factual corrections to committed files. The rest are scope decisions.
