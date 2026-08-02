# Project Ideas — High-Schooler Computational Biophysics on an SGE GPU Cluster

A curated set of **30 fully-fleshed research projects** designed to make maximal use of an
SGE/OGS cluster like the BU SCC (job arrays via `qsub -t`, GPU nodes requested by compute
capability, checkpointing, `-pe omp` threading).

**Two series.**
- **Ideas 01–15** — the *physics* series: MD, FEP, docking, scanning, coarse-graining
  (cluster = the compute engine).
- **Ideas 16–30** — the *biophysics × data-science* series: ML, statistics, information
  theory, and uncertainty applied to stability/FEP/MD (cluster = the data-generator and the
  sweep engine). These are the ones that make a project stand out among peers — every one
  pairs a physical question with a data-science method that most students won't have touched
  (conformal prediction, active learning, Markov state models, physics-informed ML,
  protein-language models, bootstrap audits, information-theoretic sampling).

Each idea is written so that the **cluster is load-bearing, not a luxury**: the project's
science would be impossible (or take a year on a laptop) without the embarrassingly-parallel
fan-out that SGE arrays provide. That is the whole point of a cluster project — you split one
big computation into thousands of independent, resumable jobs.

**How to read each idea file.** Every file follows the same template:

1. **Research question** — the thing you want to answer (competition-worthy, publishable).
2. **Why this needs an SGE cluster** — the fan-out geometry (array size = × tasks) and which
   cluster features it uses (arrays / GPU / threading / checkpointing).
3. **Data & tools** — free datasets, software, and licenses.
4. **Skill prerequisites** — what to learn first; the difficulty ladder entry.
5. **Cluster budget** — rough GPU-hours / task-hours so you can size the run.
6. **Milestones** — a staged plan from "day one" to "presentable result".
7. **Deliverables** — figures/tables a judge can hold.
8. **Pitfalls** — the traps (learned the hard way in the parent SOD1-FEP project).

---

## The SGE-usage matrix

| # | Idea | Fan-out geometry | GPUs | Arrays | Threads | Checkpoint | Reuses existing stack |
|---|------|------------------|------|--------|---------|-----------|----------------------|
| 01 | In-silico saturation mutagenesis scan | ~3,000–5,000 single mutations | No | ★★★ | ★ | – | FoldX/Rosetta |
| 02 | Consensus force-field FEP benchmark | variants × 4 FFs × 3 reps | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |
| 03 | Batched ColabFold proteome scan | hundreds of proteins | ★★★ | ★★★ | – | ★★ | ColabFold/AlphaFold |
| 04 | Virtual drug-repurposing screen | thousands of ligands | ★★ | ★★★ | ★ | ★★ | Vina/GNINA |
| 05 | Deep-mutational-scan validation | ~1,000 mutations | – | ★★★ | ★ | – | FoldX/Rosetta |
| 06 | Binding free-energy consensus | ligands × reps × FFs | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |
| 07 | Peptide-library interface scan | thousands of peptides | ★★ | ★★★ | ★ | ★★ | Vina + GROMACS |
| 08 | Coarse-grained self-assembly sweep | compositions × temperatures × seeds | ★★ | ★★★ | ★★ | ★★★ | GROMACS + MARTINI |
| 09 | Kinetic stability from unfolding replicates | variants × replicates | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |
| 10 | Dimerization free-energy panel | variants × legs × windows | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |
| 11 | Thermostability redesign | candidate mutations × FFs | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |
| 12 | Cryptic-pocket detection ensemble | replicas × pocket frames | ★★★ | ★★★ | ★ | ★★★ | GROMACS + fpocket |
| 13 | Host–pathogen docking scan | pathogen × host proteins | – | ★★★ | ★ | ★ | HADDOCK/Vina |
| 14 | ML surrogate on cluster data | hyperparameter grid | ★★ | ★★★ | ★ | ★★ | scikit-learn/GNN |
| 15 | Soft-core protocol benchmark | mutations × soft-core variants | ★★★ | ★★★ | ★★ | ★★★ | **GROMACS + pmx (this repo)** |

### Ideas 16–30 — biophysics × data science

| # | Idea | Fan-out geometry | GPUs | Arrays | Threads | Checkpoint | Method family |
|---|------|------------------|------|--------|---------|-----------|---------------|
| 16 | Conformal prediction for ΔΔG | calibration scans × predictors | – | ★★★ | ★ | – | Uncertainty quantification |
| 17 | Active learning for the gate | rounds × FEP windows | ★★★ | ★★★ | ★ | ★★★ | Experimental design |
| 18 | Cycle-consistency physics-informed ML | label scan + cycle/triple enumeration | – | ★★★ | ★ | – | Physics-informed ML |
| 19 | Markov state models of folding | trajectories (60+/variant) | ★★★ | ★★★ | ★ | ★★★ | Kinetic networks |
| 20 | Learned collective variables (VAMPnets) | trajectory ensemble + training sweep | ★★★ | ★★★ | ★ | ★★★ | Representation learning |
| 21 | Epistasis detection from DMS | double-mutant pairs (~N²) | – | ★★★ | ★ | – | Interaction modeling |
| 22 | Protein-language-model embeddings | sequence scoring + model sweep | ★★★ | ★★★ | ★ | ★★ | Representation learning |
| 23 | Predictability ceiling (noise floor) | bootstrap / MC resamples | – | ★★★ | ★ | – | Information theory |
| 24 | Transfer learning across proteins | source scans + train/test grid | – | ★★★ | ★ | ★★ | Transfer learning |
| 25 | Bootstrap audit of MBAR errors | window refits (thousands) | ★★★ | ★★★ | ★ | ★★★ | Uncertainty audit |
| 26 | Predicting FEP window failure | probe runs + label array | ★★★ | ★★★ | ★ | ★★★ | Applied reliability ML |
| 27 | Co-evolution coupling vs FEP | FEP panel + DCA pairs | ★★★ | ★★★ | ★ | ★★★ | Multi-signal fusion |
| 28 | Gaussian-process λ-ladder | adaptive windows × repeats | ★★★ | ★★★ | ★ | ★★★ | Active learning / GP |
| 29 | Information-theoretic frame selection | window selection experiments | ★★★ | ★★★ | ★ | ★★★ | Information theory |
| 30 | Anomaly detection on MD trajectories | trajectories + scoring passes | ★★★ | ★★★ | ★ | ★★★ | Unsupervised ML |

Legend: ★ = minor use, ★★ = significant, ★★★ = essential.

---

## Difficulty ladder

```
Easiest → hardest
───────────────────────────────────────────────────────────────
 03  ColabFold batch          (beginner: run scripts, interpret)
 07  Peptide docking scan     (beginner-intermediate)
 04  Drug-repurposing screen  (beginner-intermediate)
 01  Saturation scan          (intermediate: pipeline plumbing)
 05  DMS validation           (intermediate: stats + benchmarks)
 13  Host–pathogen docking    (intermediate)
 12  Cryptic-pocket ensemble  (intermediate)
 08  CG self-assembly sweep   (intermediate: MD know-how)
 09  Kinetic stability        (intermediate-advanced: MD + stats)
 14  ML surrogate             (intermediate-advanced: ML + data)
 11  Thermostability redesign (advanced: FEP)
 06  Binding FF consensus     (advanced: FEP + MM/PBSA)
 10  Dimerization FEP panel   (advanced: FEP, hardest physics)
 02  Consensus-FF benchmark   (advanced: FEP at scale)
 15  Soft-core benchmark      (advanced: FEP protocol research)
───────────────────────────────────────────────────────────────
```

**Series 2 (16–30) — biophysics × data science difficulty ladder:**

```
Easiest → hardest
───────────────────────────────────────────────────────────────
 16  Conformal prediction      (intermediate: stats-heavy, no MD)
 23  Predictability ceiling    (intermediate: stats / meta-analysis)
 21  Epistasis from DMS        (intermediate: fits + structure)
 22  pLM embeddings            (intermediate-advanced: ML + GPU)
 26  FEP failure prediction    (intermediate-advanced: applied ML)
 27  Co-evolution vs FEP       (intermediate-advanced: multi-signal)
 30  Anomaly detection         (intermediate-advanced: MD + ML)
 24  Transfer learning         (advanced: ML methodology)
 17  Active learning gate      (advanced: experimental design)
 25  Bootstrap MBAR audit      (advanced: simulation statistics)
 18  Cycle-consistency ML      (advanced: physics-informed ML)
 19  Markov state models       (advanced: kinetics + stats)
 28  GP λ-ladder               (advanced: active learning + FEP)
 29  Frame selection           (advanced: info theory + FEP)
 20  VAMPnets                  (advanced: deep learning + MD)
───────────────────────────────────────────────────────────────
```

**Recommended path for a first project:** start with **03** or **04** (fast win, single
`qsub -t` array, teaches the cluster), then do **01** or **05** (adds statistics and a
validation story), and only then attempt an FEP project (**02/06/09/10/11/15**) once you can
read `dhdl.xvg` output and debug a failed `mdrun`.

**For the data-science track:** **23** or **16** are the friendliest no-MD entry points
(pure stats with biophysics meaning), then **22** or **26**, then **17/18/19/28** once you
can run the parent repo's FEP and read MBAR output.

---

## Rough cluster budgets

These assume an 8-GPU node pool and ~12 h walltime per task (the BU SCC buy-in nodes).

| Project | Total jobs | Per-job time | Wall-clock on 8 GPUs |
|---------|-----------|--------------|----------------------|
| 01 saturation scan (FoldX) | 3,800 | ~5 min CPU | ~2 days (CPU, threaded) |
| 02 consensus FEP (8 gate variants) | 8×4×3×2 = 192 windows | ~3–6 h GPU | ~1 week |
| 03 ColabFold (200 proteins) | 200 | ~15–45 min GPU | ~1–2 days |
| 04 drug screen (3,000 ligands) | 3,000 | ~10 min CPU | ~3 days (CPU) |
| 05 DMS validation (1,000 muts) | 1,000 | ~5 min CPU | ~1–2 days |
| 06 binding FF consensus (5 ligands) | 5×3×4×3 = 180 windows | ~3–6 h GPU | ~1 week |
| 07 peptide scan (2,000 peptides) | 2,000 | ~15 min CPU | ~2–3 days |
| 08 CG sweep (64 systems × 5 seeds) | 320 | ~6 h GPU | ~2–3 days |
| 09 kinetic stability (10 variants × 5 reps) | 50 | ~12 h GPU | ~3–4 days |
| 10 dimer FEP (10 variants × 2 legs × 18 × 3) | 1,080 windows | ~3–6 h GPU | ~2–3 weeks |
| 11 thermostability redesign (20 candidates) | 20×3×2 = 120 windows | ~3–6 h GPU | ~1 week |
| 12 cryptic pockets (32 replicas) | 32 | ~12 h GPU | ~2 days |
| 13 host–pathogen docking (500 pairs) | 500 | ~1–2 h CPU | ~2–3 days |
| 14 ML surrogate (500 configs) | 500 | ~20 min GPU | ~1–2 days |
| 15 soft-core benchmark (12 windows × 2 sc) | 24 | ~3–6 h GPU | ~1–2 days |

Series 2 budgets (16–30; several are ML/stats-lean but need a *data-generating* array first):

| Project | Total jobs | Per-job time | Wall-clock on 8 GPUs |
|---------|-----------|--------------|----------------------|
| 16 conformal calibration (FoldX scan) | ~2,000 | ~5 min CPU | ~1–2 days |
| 17 active learning (10 rounds × 2 variants) | ~2,160 FEP windows | ~3–6 h GPU | ~2–4 weeks |
| 18 cycle-consistency ML (scan + training) | ~3,000 labels + ~200 configs | ~5 min CPU / ~20 min GPU | ~2–3 days |
| 19 MSM (4 systems × 60 × 150 ns) | 240 GPU runs | ~6 h GPU | ~1–2 weeks |
| 20 VAMPnets (shared trajectories + sweep) | 240 GPU + ~100 configs | ~6 h GPU / ~30 min GPU | ~1–2 weeks |
| 21 epistasis (double-mutant pairs) | ~5,000–20,000 | ~5–10 min CPU | ~2–3 days |
| 22 pLM embeddings (family scan + models) | ~10⁴ seqs + ~100 configs | min GPU | ~2–3 days |
| 23 predictability ceiling (resamples) | ~5,000 | sec–min CPU | ~1 day |
| 24 transfer learning (sources + grid) | ~10⁴ labels + ~600 runs | ~5 min CPU / ~20 min GPU | ~2–4 days |
| 25 MBAR bootstrap audit (5-rep subset) | ~1,440 FEP windows + ~10⁴ refits | ~3–6 h GPU | ~1–2 weeks |
| 26 FEP failure prediction (probe arrays) | ~1,000 (500 real + 500 probe) | ~10 min–6 h | ~1–2 weeks |
| 27 co-evolution vs FEP (20-variant panel) | ~2,200 FEP windows + DCA | ~3–6 h GPU | ~2–3 weeks |
| 28 GP λ-ladder (6–8 variants) | ~1,000–1,500 windows | ~3–6 h GPU | ~1–2 weeks |
| 29 frame selection (long windows) | ~72 windows + ~10⁴ refits | ~6 h GPU | ~1 week |
| 30 anomaly detection (4 systems × 30–60) | ~120–240 GPU runs | ~6 h GPU | ~2–4 weeks |

---

## How to pick (selection guide)

- **Only have 3 months?** → 03, 04, 01, or 05. All have a clean "dataset → number → plot"
  arc and need only CPU or modest GPU.
- **Have a semester + a mentor with MD experience?** → 08, 09, 12, or 06.
- **Want the strongest "real science" story (ISEF / Regeneron STS / JSS)?** → 02, 10, 11, or
  15 — these produce an actual thermodynamic result with error bars, the gold standard for a
  computational-physics project. They are also the hardest.
- **Want to reuse what this repo already has?** → 02, 06, 09, 10, 11, 15 all sit directly on
  the GROMACS + pmx engine, `scripts/submit_array.sh` array pattern, and the
  provenance-stamped NPZ + MBAR analysis already in `src/`.
- **Want to stand out with data science?** → the 16–30 series is built for this. The best
  "new idea" combinations that still reuse the repo's data:
  - **25 + 02** — bootstrap-audit the consensus-FF errors (rigorous methods + your own data).
  - **16 + 01** — conformal intervals around a FoldX saturation scan.
  - **19 + 20** — MSM then VAMPnets on the same Stage-4 trajectory ensemble (two analyses, one
    data-generating array).
  - **28 + 02** — GP-guided λ placement evaluated inside the consensus-FF benchmark.
  - **26 + 10** — predict failures before paying for the dimer-FEP array.
  - **22 + 05** — pLM embeddings vs physics on a DMS (the "pLM vs physics" debate, tested).
  - **17 + 02** — active learning to choose *which* controls get FEP'd first.

---

## Cluster mechanics cheat-sheet (from the parent project's hard-won lessons)

- **One array task = one (variant/ligand/mutation, leg, window, replicate).** Do not merge
  dimensions into one long job — resumability and failure isolation depend on it.
- Request GPUs by compute capability, not model: `-l gpus=1 -l gpu_c=7.0`. The scheduler sets
  `CUDA_VISIBLE_DEVICES`; never hardcode a device ID.
- Keep `-l h_rt` ≤ 12 h to stay eligible for buy-in nodes.
- `source scripts/scc_env.sh` is the ONE way to get a working shell (module chain + conda env).
  Batch scripts must source it too or every array task dies identically.
- **Checkpoint everything.** A window that dies at hour 11 must resume, not restart.
- Do not persist raw trajectories when you only need estimates — disk fills fast.
- Test with a 3-task `qsub -t 1-3` array before launching the full 1,000+.

---

## SGE quick reference

```bash
# submit an array of N tasks
qsub -t 1-N -v MY_VAR=x scripts/your_job.sh

# decode the task id inside the job
SGE_TASK_ID   # 1..N
SGE_TASK_FIRST, SGE_TASK_LAST, SGE_TASK_STEP

# status / hold / delete
qstat -u $USER          # or qstat -u $USER -j <jobid>
qhold <jobid>           # pause
qdel <jobid>            # kill

# interactive GPU node
qrsh -l gpus=1 -l gpu_c=7.0 -pe omp 8

# project charge
qsub -P <project> ...
```
