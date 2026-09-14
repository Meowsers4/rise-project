# SOD1 Variant Stability Pipeline

A completed GPU-cluster validation campaign for protein-stability changes (ΔΔG of
folding) in clinically observed **SOD1** variants using **alchemical free-energy
(FEP/TI)** calculations. The failed gate stopped the planned extension to
uncharacterized ALS variants; the project is now a methods-and-limitations study.

> **Status (2026-09-14):** the pre-registered validation gate failed (Pearson
> *r* = 0.326; RMSE = 2.123 kcal/mol). By explicit user decision, the project
> retains apo-2SH, ends the current GPU campaign, retires charge-changing coverage
> and uncharacterized-variant triage, and proceeds as a methods-and-limitations
> deliverable. No new GPU submission is authorized. Read §2.3 and §11 before
> writing any abstract, poster, or paper claim.

---

## 0. Agent orientation — read this first

If you are an agent resuming work on this project:

1. **The method is settled.** GROMACS + pmx. Do not re-open the framework decision.
   Do not port to OpenMM/Perses or AMBER. See [§2.4](#24-framework-gromacs--pmx).
2. **The novelty story changed.** This project is no longer "FEP reveals SOD1
   variants destabilize" — that was published in 2021 by Wells et al. using the
   *same toolchain*. Read [§11](#11-prior-art-and-novelty-position) for what is
   dead, what survives, and what the defensible claim now is.
3. **The 2026-09-14 scope pivot is binding.** See [§2.3](#23-claim-contract-and-approved-pivot):
   C1 and C3 are retired, C4 is deferred and not claimed, and C2 survives only as
   a methods-and-limitations result.
4. **No new GPU work is authorized.** The remaining work is analysis, archival
   verification, and preparation of the methods-and-limitations deliverable.
5. **Changelog is [§12](#12-changelog).** Append to it; don't rewrite history.

---

## 1. What this project is (and is not)

**Original goal:** produce a validated ΔΔG map across the SOD1 variant landscape.
The validation gate did not support that use. The approved current goal is to
report the failed gate, sampling/convergence diagnostics, benchmark chain of
custody, reference-state sensitivity, and the negative G93A disulfide diagnostic
as a reproducible methods-and-limitations study.

This is now a **methods/target limitations** project:
- *Target:* SOD1 (Cu/Zn superoxide dismutase, UniProt **P00441**), a 153-residue
  soluble homodimer implicated in familial ALS.
- *Method question:* what do the failed predictive gate, convergence diagnostics,
  and reference-state mismatch show about apo-2SH SOD1 stability FEP?

**In scope now:** reproducible analysis of completed control calculations,
benchmark/reference-state audit, the G93A disulfide negative test, archival
verification, and a methods-and-limitations deliverable.

**Out of scope for the current project:** new GPU FEP/MD, charge-changing support,
VUS prediction/classification, concordance claims, wet-lab work, and holo/metal-
bound or apo-SS campaigns. A later arm requires a new explicit scope decision and
preregistration.

---

## 2. Core design decisions (do not silently override)

These decisions propagate through the entire pipeline. An agent modifying this
project must treat changing any of them as a scope change, not a refactor.

### 2.1 Apo-first
Simulate the **apo, disulfide-reduced (2SH)** form in v1.
- **Why:** avoids metalloprotein force-field parameterization (Cu/Zn coordination),
  AND lands on the aggregation-prone, disease-relevant species.
- **Novelty caveat (new):** this decision is *scientifically correct but not novel*.
  Both Wells 2021 and Hsueh 2022 simulated apo-2SH. Do **not** frame apo-first as a
  contribution. Its value here is that it is the state where Wells reported
  **convergence failure** — which is the opening. See [§11.3](#113-what-survives-as-novelty).
- Holo/metal-bound is a **separate parameterized arm**, not a config flag.

### 2.2 Control validation is a GATE, not a final step
Before spending cluster time on uncharacterized variants, reproduce the **known
experimental ΔΔG** of the characterized controls (A4V, G93A, G37R, …).
- If the control correlation is poor, extending to novel variants produces
  confident nonsense. Treat this as **go/no-go**.
- **Pre-registered thresholds (new — set these before looking at results):**
  - Pearson r ≥ 0.70 against the Kumar 2017 control set, **and**
  - Spearman ρ reported alongside (rank agreement matters more than absolute
    scale for triage), **and**
  - median |cycle-closure hysteresis| below a stated kcal/mol bound (fix the bound
    in `config/pipeline.yaml` before the first gate evaluation).
- Reference point: Wells 2021 reported correlation **0.81** with experiment on 10
  variants. Matching or beating that on a 54-variant panel is the quotable result.
- The control-reproduction result is itself a presentable outcome — and under the
  revised framing, it is a *primary* deliverable, not a preamble.

### 2.3 Claim contract and approved pivot

The original four claims remain recorded for chain of custody. After the failed
validation gate and reference-state audit, the user approved this scope pivot on
2026-09-14. Downstream sections describing charge-changing or VUS stages are
historical design notes; this decision governs wherever they conflict.

| # | Original claim | Status after sign-off |
|---|---|---|
| **C1** | Extension to **charge-changing** variants that Wells 2021 excluded | **Retired.** No charge-changing implementation or GPU work in the current project. |
| **C2** | **Convergence of apo-2SH** where Wells reported non-convergence | **Retained only as a methods-and-limitations result.** Report the diagnostics, failed predictive gate, and reference-state caveats; do not claim validated predictive convergence. |
| **C3** | **Prospective, DMS-independent VUS triage** | **Retired.** No VUS FEP predictions or classifications in the current project. |
| **C4** | **Concordance/discordance analysis** vs DMS + FoldX/Rosetta/ML | **Deferred and not claimed.** It has no current-project FEP input once C3 is retired. |

The current load-bearing output is the reproducible benchmark audit,
convergence/uncertainty diagnostics, reference-state sensitivity analysis, and
negative G93A disulfide test.

**Do not write "first," "novel FEP protocol," or "we show SOD1 variants are
destabilizing."** All three are false. See [§11.2](#112-what-is-dead).

### 2.4 Framework: GROMACS + pmx
**Resolved.** GROMACS + pmx, equilibrium λ-windows with **MBAR** estimation.

- **This is the same toolchain as Wells 2021.** That is a deliberate choice, not an
  oversight: it makes the comparison apples-to-apples rather than confounded by
  engine differences.
- **But it forecloses one novelty avenue.** An earlier draft of this project
  considered OpenMM+Perses, where "first rigorous folding-stability validation of
  Perses" would have been a real methods contribution (Perses is documented for
  relative *binding* free energy and protein:protein binding, and is conspicuously
  absent from the folding-stability benchmark literature). **With pmx, that claim is
  unavailable** — pmx/GROMACS is the most extensively benchmarked stability-FEP
  stack in existence.
- **Estimator difference is still real but thin.** Wells used pmx in
  **non-equilibrium fast-growth** mode (Crooks/Jarzynski work analysis, ~50 ps
  transitions seeded from long equilibrium trajectories). This project uses
  **equilibrium windowed sampling + MBAR**. That is a genuinely different estimator
  and supports C2, but it is not by itself a headline.

---

## 3. Pipeline overview

```
Stage 0  Build variant panel        (CPU, data)     controls + uncharacterized
Stage 1  Prepare structure          (CPU)           apo, protonate, solvate
Stage 2  Cheap prescreen            (CPU)           FoldX / Rosetta, ALL variants
Stage 3  Alchemical FEP  ★          (GPU, cluster)  ΔΔG — the workhorse
Stage 4  Unbiased MD                (GPU)           mechanism, subset only
Stage 5  Validate + concordance     (CPU)           controls gate, then triage
```

★ Stage 3 is the parallel cluster stage. Its job count fans out as:

```
jobs  =  variants × 2 legs × λ-windows × replicates
      ≈  40 × 2 × 18 × 5  ≈  7,200 independent GPU jobs
```

Every job is embarrassingly parallel (no inter-job communication), so throughput
scales linearly with GPU count.

---

## 4. Stages in detail

### Stage 0 — Build and stratify the variant panel
- Sources: **ClinVar**, **ALSoD** (ALS-specific), cross-ref **gnomAD** for frequency.
- **Numbering caution:** mature SOD1 has the initiator Met removed → 153 residues.
  Literature variant names (A4V, G93A) use mature numbering. Get the offset right
  once, centrally, or every downstream mapping breaks.
- Stratify into four buckets:
  1. **Positive controls** — experimentally measured ΔΔG / ΔTm (non-negotiable).
  2. **Negative controls** — common/benign, should not destabilize.
  3. **Pathogenic-but-uncharacterized** — clinical label, no structural work.
  4. **VUS** — variants of uncertain significance (the payoff targets).
- **New required columns (audit consequence):**
  - `charge_change` — bool. Whether WT→mut alters net charge. **Wells 2021 excluded
    all of these.** Flag them explicitly; they carry claim C1.
  - `wells2021` — bool. Whether this variant appears in the Wells 10
    (A4V, A4S, A4T, L38V, G41S, G93A, G93S, I113T, V148G, V148I). These are the
    head-to-head comparison set.
  - `axakova_class` — the Axakova 2025 DMS call (abundance and function scores)
    where available. **Blocks the VUS claim until populated** — see
    [§9](#9-open-decisions).
- Output: `data/variants.csv` with columns
  `variant,mature_pos,wt_aa,mut_aa,bucket,exp_ddg,exp_source,clinvar_id,charge_change,wells2021,axakova_class`.

### Stage 1 — Prepare structure
- Start from a high-resolution human SOD1 crystal structure (PDB).
- **Apo** (see §2.1). Strip metals; model disulfide-reduced Cys57/Cys146.
- Monomer vs dimer: **dimer** for interface-adjacent variants, **monomer** for
  buried-core variants — decided per variant, recorded in `variants.csv`.
- Standard prep: fix missing atoms/loops, assign protonation (PROPKA or H++),
  solvate, neutralize, add ions.
- Deliver as a **parameterized script**: `variant_id -> prepared_system`. It runs
  dozens of times; no manual steps.

### Stage 2 — Cheap empirical prescreen (CPU, runs on everything)
- **FoldX BuildModel** and/or **Rosetta cartesian_ddg** on ALL variants.
- **Reframed:** these are no longer just a prioritization step. Under claim C4 they
  are **comparators** — the analysis of *where cheap predictors fail and why*, with
  FEP and DMS as reference, is a deliverable in its own right (and the natural
  MLSB/PSB angle).
- Add ML comparators to the same table: DDGun, DynaMut2, mCSM, RaSP,
  ThermoNet/Stability Oracle where runnable, plus AlphaMissense as a
  pathogenicity (not stability) axis.
- Validate these against controls too — if FoldX can't reproduce known cases,
  you learn it cheaply.

### Stage 3 — Alchemical ΔΔG (historical design; GPU campaign ended) ★

The items below document the original plan and completed implementation. The
2026-09-14 decision retired C1 and ended GPU work; they are not instructions to
submit further arrays.
- Thermodynamic cycle: run the alchemical mutation (X→Y) in the **folded protein**
  and in an **unfolded reference** (solvated capped tripeptide). Difference of the
  two legs = ΔΔG_folding.
- Each leg → **12–24 λ-windows**, softcore potentials for appearing/disappearing
  atoms, **≥3 (ideally 5) replicates** per window with distinct seeds.
- Replicate discipline is where people cut corners; error bars are mandatory here.
- **Framework:** GROMACS + pmx (resolved — see §2.4). Hybrid topology generation via
  `pmx mutate` / `pmx gentop`.
- Free-energy estimator: **MBAR** (pymbar). Check convergence (forward/backward
  hysteresis, cycle closure).
- **Charge-changing mutations (claim C1) — the technically hard part.** Net-charge
  changes under PME introduce finite-size and net-charge artifacts. This is exactly
  why Wells 2021 restricted to charge-conserving variants. Options to evaluate and
  document:
  - co-alchemical counterion (transform a distant ion in the opposite direction so
    total charge is conserved), or
  - Rocklin-type analytical finite-size / net-charge corrections applied post hoc.
  Whichever is chosen, it must be **validated on charge-changing controls that have
  experimental ΔΔG** before being trusted on charge-changing VUS. Treat this as a
  sub-gate inside the main gate.
- **Convergence logging (claim C2) is mandatory, not optional.** Per-window overlap
  matrices, forward/backward hysteresis, replicate spread, and time-series of the
  running estimate must be persisted for every job. C2 cannot be claimed
  retroactively from results that weren't instrumented.

### Stage 4 — Unbiased MD for mechanism (retired from current scope)
- FEP gives a number; MD gives the story. Run longer replicate plain-MD on the
  subset FEP flags as interesting — prioritize **discordant** variants from Stage 5
  (where FEP and DMS disagree), since those are where mechanism is most informative.
- Analyze: per-residue **RMSF**, contact-map / salt-bridge loss, **SASA** change,
  dimer-interface integrity, local-unfolding proxies.
- **Caveat:** do not expect full unfolding in accessible timescales — read proxies
  and equilibrium shifts, not a complete denaturation event.

### Stage 5 — Validation, triage, and concordance (closed after gate failure)

The gate was evaluated and failed. The frozen rules below remain as design
history; the approved scope does not proceed to charge or VUS stages.
- **GATE:** correlate computed FEP ΔΔG vs experimental controls against the
  pre-registered thresholds in §2.2. Go/no-go.
- **Head-to-head:** report the Wells-10 subset separately. Same target, same engine,
  different estimator, more sampling — this is the cleanest possible comparison and
  directly supports C2.
- **Charge-changing sub-gate:** report charge-changing controls separately from
  charge-conserving ones. If the charge handling doesn't validate, C1 is
  unsupported and charge-changing VUS results must be withheld.
- **Concordance table (claim C4) — the centerpiece.** For each VUS, one row:
  FEP ΔΔG ± error │ Axakova DMS abundance │ Axakova DMS function │ FoldX │ Rosetta │
  ML predictors │ AlphaMissense │ ClinVar status.
  - Lead with **agreement** as multi-line evidence (ACMG-style framing: FEP is
    physics-based and training-data-free, therefore methodologically independent of
    both DMS and of ML predictors trained on ΔΔG datasets).
  - Then lead *harder* with **disagreement**. Discordant variants are the most
    defensible new content in the project — e.g. stability-neutral but
    function-disrupting variants, which Axakova's own abundance-vs-function split
    already exposes. Do not bury these.
- Output: `results/ddg_map.csv`, `results/concordance.csv`, validation plots.

---

## 5. Cluster orchestration

- **Workflow manager:** Snakemake or Nextflow — the DAG
  (prep → prescreen → FEP windows → analysis) must be reproducible and restartable.
- **Scheduler:** SGE/OGS (BU SCC) **job arrays** (`qsub -t`) — one array task per
  `(variant, leg, window, replicate)`. GPUs requested by compute capability
  (`-l gpus=1 -l gpu_c=7.0`), never by model; the scheduler assigns the device via
  `CUDA_VISIBLE_DEVICES`, so the FEP code must not hardcode it.
- **Checkpoint everything** — windows die; never restart from zero.
- **Trajectory retention policy up front** — full FEP trajectories are large.
  Keep free-energy estimates + analysis outputs; downsample or discard raw frames.
  **Exception:** retain enough per-window data to reconstruct convergence
  diagnostics, since claim C2 depends on them.

---

## 6. Proposed repository layout

```
sod1-fep/
├── README.md                  # this file
├── env/
│   ├── environment.yml        # conda: gromacs, pmx, mdtraj, pymbar, etc.
│   └── modules.md             # cluster module-load notes
├── data/
│   ├── variants.csv           # Stage 0 output (the source of truth)
│   ├── axakova_dms.csv        # Axakova 2025 DMS scores, for concordance
│   └── structures/            # prepared PDBs, per variant
├── config/
│   └── pipeline.yaml          # global params incl. cluster block (SGE) + gate thresholds
├── workflow/
│   └── Snakefile              # or main.nf
├── src/
│   ├── panel/                 # Stage 0: pull + stratify variants
│   ├── prep/                  # Stage 1: variant_id -> prepared_system
│   ├── prescreen/             # Stage 2: FoldX / Rosetta / ML comparator wrappers
│   ├── fep/                   # Stage 3: pmx setup, run, MBAR analysis
│   ├── md/                    # Stage 4: unbiased MD + analysis
│   └── analysis/              # Stage 5: validation, concordance, plots
├── scripts/
│   └── submit_array.sh        # SGE array submission helper (qsub -t)
├── results/
│   ├── ddg_map.csv
│   ├── concordance.csv
│   ├── convergence/           # per-window diagnostics — claim C2 depends on these
│   └── figures/
└── tests/
    └── test_numbering.py      # guard the mature-numbering offset
```

---

## 7. Suggested tech stack

| Concern            | Choice                                              |
|--------------------|-----------------------------------------------------|
| Language           | Python 3.11+                                         |
| MD / FEP engine    | **GROMACS + pmx** (resolved)                         |
| Empirical ΔΔG      | FoldX, Rosetta cartesian_ddg                         |
| ML comparators     | DDGun, DynaMut2, mCSM, RaSP, AlphaMissense           |
| Free-energy math   | pymbar (MBAR)                                        |
| Trajectory analysis| MDTraj / MDAnalysis                                  |
| Workflow           | Snakemake or Nextflow                                |
| Scheduler          | SGE/OGS (BU SCC), job arrays (`qsub -t`)              |
| Structure prep     | PDBFixer, PROPKA/H++                                 |

---

## 8. Milestones

1. **M0** — variant panel built and stratified; numbering test passes;
   `charge_change` and `wells2021` columns populated.
2. **M0.5** — Axakova 2025 DMS scores ingested; VUS intersection computed.
   **Blocks claim C3.**
3. **M1** — structure-prep script runs end-to-end on one control variant.
4. **M2** — empirical prescreen + ML comparators over all variants; ranked table.
5. **M3** — FEP runs for the control set; **validation gate** evaluated against
   pre-registered thresholds. Charge-changing sub-gate evaluated separately.
   **Completed: gate failed; no further compute.** The charge-changing sub-gate was
   retired before implementation.
6. **M3.5** — Wells-10 head-to-head reported; convergence diagnostics for apo-2SH
   assembled. **Current methods-and-limitations deliverable.**
7. ~~**M4** — FEP extended to uncharacterized variants; ΔΔG map produced.~~
   **Retired with C3.**
8. ~~**M5** — concordance table built; discordant variants identified.~~
   **Deferred; C4 is not claimed.**
9. ~~**M6** — mechanism MD on discordant subset.~~ **Retired.** Final report and
   figures proceed using the completed control, audit, and diagnostic evidence.

---

## 9. Open decisions

### Resolved
- ~~FEP framework~~ → **GROMACS + pmx** (see §2.4).
- ~~λ-window count~~ → **20** (18 until `927b302`, which refined the endpoint; the changelog
  entry for 2026-08-07 records the superseded value of 18 and is left as history).
  ~~Replicate count~~ → **3 in `config/pipeline.yaml`.** The former 3 → 5 plan was retired
  when the GPU campaign ended on 2026-09-14. **Do not change config or array bounds.**

### Retired or deferred by the 2026-09-14 approved pivot
- **Charge-changing handling:** not selected or implemented. C1 is retired.
- **Axakova intersection:** exactly how many of the target VUS does Axakova 2025
  already resolve? The panel currently assumes ~38 VUS; Axakova reports 156 SOD1
  missense variants in ClinVar of which ~26% (~41) are VUS, and provides evidence
  for 41% of previously-reported VUS. **C3 is retired and C4 is deferred, so this
  is not current-project work.**
- **Gate threshold values:** were frozen before evaluation; the gate failed and remains
  the gate-of-record.

### Still open for closeout
- **SCC archive verification:** confirm that live `results/fep/G93A` is the 2SH baseline
  and that the SS diagnostic has a second SCC copy.
- **LiveCoMS inquiry:** fill author identity, expertise, and stewardship fields; approve
  the license/update plan; then send manually.
- **Independent reproducibility check:** a second person should reproduce Table 3 and
  Figure 4 from the manifest and frozen inputs.

### Closed — do not reopen
- **Alternative target fallback (TTR):** **rejected.** TTR FEP/TI stability is
  already published (His88 mutations by TI, PMID 35484710, "excellent agreement"
  with experiment), a 2025 integrative computational mutational-landscape paper
  exists (npj Syst Biol Appl, doi 10.1038/s41540-025-00582-2), and TTR adds
  homotetramer cost. It removes the metal problem but offers **no cleaner novelty
  story than SOD1**. Only revisit if SOD1 validation fails outright.

---

## 10. Pivot triggers

Explicit conditions under which the framing must change. An agent hitting one of
these should stop and surface it rather than proceeding.

| Trigger | Response |
|---|---|
| Gate Pearson r < 0.60, or large hysteresis | **Triggered and adopted 2026-09-14.** Do not extend to VUS; proceed with the methods/sampling-limits deliverable. |
| Axakova already confidently classifies most target VUS | Drop "reclassification" framing. Go all-in on C4: physics-vs-DMS-vs-ML concordance and discordance. |
| apo-2SH won't converge within GPU budget | Make convergence *the* paper. Restrict biological claims to folded/holo states. C2 becomes the whole contribution, stated as a limit. |
| Charge-changing sub-gate fails | Withhold all charge-changing VUS results. C1 is unsupported; report the failure honestly — it corroborates Wells' reason for excluding them. |

---

## 11. Prior art and novelty position

*Added 2026-08-07 following a literature audit. This section exists so that no one
— human or agent — reconstructs a novelty claim that was already checked and found
to be occupied.*

### 11.1 The direct precedent

**Wells NGM, Tillinghast GA, O'Neil AL, Smith CA. "Free energy calculations of
ALS-causing SOD1 mutants reveal common perturbations to stability and dynamics
along the maturation pathway." *Protein Science* 2021;30(9):1804–1817.**
doi [10.1002/pro.4132](https://doi.org/10.1002/pro.4132) · PMID 34076319 · PMC8376412

- **Same target, same tool (pmx), same apo-2SH framing.**
- 10 **charge-conserving** variants: A4V, A4S, A4T, L38V, G41S, G93A, G93S, I113T,
  V148G, V148I.
- **Non-equilibrium fast-growth** alchemy (Crooks/Jarzynski work analysis) — ~50 ps
  morphing transitions launched at intervals from long equilibrium trajectories,
  forward + reverse.
- States covered: apo-2SH monomer and dimer, apo-SS, Zn-only intermediate, Cu/Zn
  holo monomer, holo dimer — i.e. the full maturation pathway.
- Correlation with experiment: **0.81** (0.46 with truncated sampling).
- **Stated limitations — these are the openings:**
  1. Several residues in less-mature states **failed to converge within 250 ns**
     (undersampling of exactly the apo species this project targets).
  2. **Charge-conserving restriction** due to PME net-charge artifacts.
  3. Temperature mismatch with experiment.

> Verification note: the full text was not directly retrievable (publisher
> bot-detection). The variant roster and protocol were reconstructed from extensive
> verbatim search-index snippets and are internally consistent, but the exact
> variant table and trajectory-length description **should be checked against the
> PDF** before being cited in a manuscript. The reported 0.81 is a correlation
> coefficient; Pearson vs R² is not disambiguated in recoverable text.

### 11.2 What is dead

Do not claim any of the following. Each is anticipated:

- ❌ **"Apo-first is novel."** Anticipated by Wells 2021 *and* by Hsueh et al.,
  *Front. Mol. Biosci.* 2022, doi [10.3389/fmolb.2022.845013](https://doi.org/10.3389/fmolb.2022.845013)
  (Plotkin group — apo vs holo, 2SH vs SS, monomer vs dimer for A4V and D101N;
  found metalation adds only ~1.5 kcal/mol to dimer binding, so metal effects are
  mostly on folding, not dimerization). Apo-2SH is the *consensus* framing in this
  literature, not a differentiator.
- ❌ **"FEP reveals SOD1 variants are destabilizing."** Wells 2021.
- ❌ **"First rigorous open-source FEP stability protocol."** pmx/GROMACS is the
  most benchmarked stability-FEP stack there is: Seeliger & de Groot's barnase
  109-mutation benchmark (r≈0.86); Baby et al., *PLOS One* 2025,
  doi [10.1371/journal.pone.0335829](https://doi.org/10.1371/journal.pone.0335829)
  reports GROMACS r=0.87, AUE 0.61 kcal/mol on staphylococcal nuclease and T4
  lysozyme. Schrödinger's Protein FEP+ covers 328 mutations across 14 structures
  (Scarabelli et al., *J Mol Biol* 2022;434(2):167375,
  doi [10.1016/j.jmb.2021.167375](https://doi.org/10.1016/j.jmb.2021.167375),
  R²≈0.65, MUE≈0.95 kcal/mol).
- ❌ **Engine-validation novelty.** This *was* available under OpenMM+Perses, which
  is documented for relative binding and protein:protein binding but absent from
  the folding-stability benchmark literature (cf. Zhang et al., *JCTC* 2024 /
  bioRxiv 2023.03.07.530278, barnase:barstar, RMSE 1.61 kcal/mol on binding). It
  is **not** available under pmx. Recorded here so the tradeoff is visible if the
  framework decision is ever revisited.
- ❌ **"First SOD1 VUS reclassification."** Already attempted by DMS (Axakova 2025),
  by in-vitro aggregation propensity, and by functional/zebrafish assays
  (Bedja-Iacona et al., *IJMS* 2025, doi
  [10.3390/ijms26157414](https://doi.org/10.3390/ijms26157414), reclassifying
  p.Val120Leu as pathogenic).

### 11.3 What remains defensible after the failed gate

- **C1 — retired.** Charge-changing variants remain an interesting future methods
  problem, but this project did not implement or validate them and makes no claim.
- **C2 — methods-and-limitations only.** Wells reported apo-2SH non-convergence. This
  project can report its equilibrium-window diagnostics, failed predictive gate,
  and the limited effect of extending A4V sampling. It cannot claim that it resolved
  apo-2SH prediction or validated the method.
- **C3 — retired.** No prospective VUS FEP triage or classification will be reported.
- **C4 — deferred and not claimed.** Without VUS FEP outputs there is no preregistered
  concordance analysis to perform.
- **Benchmark/audit contribution.** The reproducible chain-of-custody audit exposes a
  heterogeneous experimental anchor, including whole-dimer G93S/G93V measurements,
  while preserving the original failed gate as the gate-of-record.
- **Negative diagnostic contribution.** Restoring the G93A disulfide shifted folded
  ΔG by only −0.017 kcal/mol, ruling out that single explanation for G93A's error
  without generalizing to other variants.

### 11.4 The orthogonal-evidence landscape

- **Axakova et al. 2025** — bioRxiv doi
  [10.1101/2025.02.25.640191](https://doi.org/10.1101/2025.02.25.640191);
  published *Am J Hum Genet* 2025;112(10):2295–2315, doi
  [10.1016/j.ajhg.2025.08.019](https://doi.org/10.1016/j.ajhg.2025.08.019).
  A **deep mutational scan**, not a computational classifier. Measured >2,000 SOD1
  substitutions (~86% of all possible missense) for **both enzymatic function and
  protein abundance**. Reports 156 SOD1 missense variants in ClinVar, ~26% VUS, and
  provides new evidence for **41% of previously-reported VUS**. The
  abundance-vs-function split is directly relevant: it separates
  stability-mediated from function-mediated effects, which is exactly the axis FEP
  speaks to.
- **Kumar et al. 2017** — *Sci Rep*, doi
  [10.1038/s41598-017-04950-9](https://doi.org/10.1038/s41598-017-04950-9). The
  experimental stability compilation used as the control anchor (~47 mutations
  across apo-monomer / holo-dimer / loopless structures; PoPMuSiC and FoldX were
  the best of 8 predictors tested; apo-monomer vs holo-dimer ΔΔG correlated
  R²≈0.70). **Not a state-homogeneous anchor:** the source audit found that its
  G93S/G93V “apo-monomer” values are whole-dimer apo DSC results at 49.4 °C, and
  all eight gate controls are disulfide-intact while this project simulates 2SH.
  See `docs/reference_state_audit.md` and `docs/reference_sensitivity_analysis.md`.
- **ML/empirical predictors already applied to SOD1:** FoldX, Rosetta
  cartesian_ddg, DDGun, DynaMut2, mCSM, SAAFEC-SEQ, MAESTRO, AlphaMissense. SOD1 is
  a textbook demonstration target — assume any given predictor has been run on it.

### 11.5 Competition framing (STS / ISEF / MIT URTC / PSB / MLSB)

Judges reward defensible rigor and a crisp "what did the test establish?" narrative
over scale. Do not say the method was validated, do not imply that it was applied
prospectively to VUS, and do not turn benchmark heterogeneity into an excuse for
the failed gate.

For the current methods-and-limitations deliverable, lead with:

1. **The pre-registered negative result.** The predictive validation gate failed,
   and the project stopped rather than extending to VUS.
2. **The failure analysis.** Separate sampling precision from reference-state and
   benchmark heterogeneity without claiming that either alone explains the error.
3. **Reproducibility.** Preserve frozen predictions, audited references, manifests,
   figures, and the negative disulfide diagnostic so another person can regenerate
   the central results.

The strongest bounded message is procedural: a preregistered gate prevented a
failed control model from being promoted into confident variant predictions, and
the follow-up diagnostics narrowed what did and did not explain that failure.

---

## 12. Changelog

### 2026-09-14 — failed-gate pivot approved
- User explicitly approved retaining apo-2SH and ending the current GPU campaign.
- C1 and C3 retired; C4 deferred and not claimed; C2 bounded to a
  methods-and-limitations result rather than predictive validation.
- The gate-of-record remains Pearson *r* = 0.326 and RMSE = 2.123 kcal/mol.
- Any apo-SS or future predictive campaign now requires a matched benchmark, new
  preregistration, and explicit approval.

### 2026-08-07 — literature audit; framework resolution; reframing
- **Framework resolved:** GROMACS + pmx (was: open between OpenMM+Perses,
  GROMACS+pmx, AMBER TI). λ-windows fixed at 18, replicates at 5.
- **Prior-art audit performed.** Wells 2021 identified as direct precedent — same
  target, same toolchain. Hsueh 2022 and Axakova 2025 identified as adjacent.
- **Project reframed** from discovery ("FEP reveals SOD1 destabilization") to a
  methods-extension + orthogonal-evidence-triage contribution. Four load-bearing
  claims formalized as C1–C4 in §2.3.
- **Perses engine-validation novelty recorded as forgone** under the pmx decision
  (§2.4, §11.2) — documented so the tradeoff is visible if revisited.
- Charge-changing variants elevated from an implementation detail to the primary
  methods claim (C1); new Stage 3 sub-gate added.
- Convergence diagnostics elevated from QC to a claim dependency (C2); retention
  policy amended in §5.
- Stage 5 rewritten around a concordance table (C4).
- Gate thresholds pre-registered (§2.2). Pivot triggers added (§10).
- TTR fallback **rejected and closed** (§9) — already modeled, no cleaner novelty.
- New `variants.csv` columns required: `charge_change`, `wells2021`,
  `axakova_class`. New milestones M0.5 and M3.5.

### Historical verification debts (not current-project blockers)
- Wells 2021 full text not directly retrievable; variant roster and protocol
  reconstructed from search-index snippets. **Verify against the PDF before citing.**
- Axakova ∩ target-VUS intersection was not computed variant-by-variant. C3 was
  subsequently retired, so this is no longer a current-project blocker.
- "No published FEP-based SOD1 VUS reclassification" is absence of evidence across
  current-literature searches. Strong, not absolute. State as such in any manuscript.

---

## 13. References for any future scope expansion

- Wells 2021 PDF — verify variant roster, λ protocol, trajectory lengths, and
  whether 0.81 is Pearson or R².
- Axakova 2025 supplementary tables — per-variant abundance and function scores,
  for `data/axakova_dms.csv`.
- Kumar 2017 supplementary — experimental ΔΔG for the control panel.
- pmx hybrid-topology setup docs and the Gapsys/de Groot protocol papers.
- Rocklin et al. finite-size / net-charge correction methodology (for C1).
- ClinVar + ALSoD current variant exports (re-pull; the 2025 literature moved).
