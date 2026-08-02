# Project Idea Summaries — 30 Research Proposals at the Intersection of Biophysics, High-Performance Computing, and Data Science

A condensed but substantive digest of all 30 projects in this directory. Each summary
captures the research question, the methodological core, why it demands an SGE cluster, what
makes it original, and the headline deliverable. Full depth lives in the individual
`idea-NN_*.md` files.

Two complementary series:

- **Ideas 01–15 — the physics series.** Molecular-dynamics, free-energy, docking, and scanning
  projects where the cluster is the *compute engine* that turns an impossible single-machine
  computation into a weekend. These reuse the parent SOD1-FEP GROMACS + pmx infrastructure
  wherever possible.
- **Ideas 16–30 — the biophysics × data-science series.** Machine learning, statistics,
  information theory, and uncertainty quantification applied to stability prediction, FEP, and
  MD trajectories. The cluster is the *data generator* (arrays of labels, trajectories,
  resamples) and the *sweep engine* (hyperparameter and experiment grids). These are the ideas
  engineered to differentiate a student from their peers: every one pairs a physical question
  with a data-science method that is rarely seen in high-school research.

---

# Series I — The Physics Series (01–15)

## Idea 01 — In-Silico Saturation Mutagenesis Scan
*Predicting protein evolvability from first principles.*

- **Question:** Map the complete stability landscape of a protein — every single-point
  mutation at every position (~19 × N ≈ 3,000–5,000 for a typical protein) — and relate it to
  evolutionary conservation.
- **Method:** An array task per mutation, each running FoldX `BuildModel` (averaged over 3–5
  runs for noise control), assembled into a position × amino-acid ΔΔG heatmap.
- **Cluster role:** ~2,900–4,000 independent CPU tasks; the 2-D landscape *only exists*
  because thousands of independent points were computed. This is the canonical
  embarrassingly-parallel workload.
- **Originality:** Most "stability" studies look at a handful of variants; a *complete*
  landscape permits the genuinely interesting claim about which regions evolution can explore
  freely — a mutability/evolvability map rather than a variant list.
- **Deliverable:** Full stability heatmap + conservation correlation + ranked intolerant and
  tolerant positions with structural explanations.

---

## Idea 02 — Consensus Force-Field FEP Benchmark
*Does averaging free energies over force fields beat any single one?*

- **Question:** Reproduce the Gapsys 2016 / Wells 2021 claim — that averaging alchemical ΔΔG
  over ≥2 force fields reduces error to ~1 kcal/mol — on your own system, using the SOD1 gate
  subset that already has experimental labels.
- **Method:** Run the parent repo's GROMACS + pmx engine under all four force fields shipped in
  `pmx/data/mutff45` (`amber99sbmut`, `amber99sb-star-ildn-mut`, `charmm22starmut`,
  `oplsaamut`); combine via mean/median and compare against experimental ΔΔG.
- **Cluster role:** 8 variants × 4 FFs × 2 legs × 18 windows × 3 reps = 3,456 GPU windows,
  each an independent 3–6 h job. Weeks of single-GPU compute becomes ~2 weeks on 8 GPUs.
- **Originality:** It is a *controlled, pre-registered methods experiment* — the only variable
  between arms is the force field — and it directly informs a config decision in the parent
  repo (which FF to trust).
- **Deliverable:** Per-FF and consensus ΔΔG tables with a bar chart showing the consensus
  landing below every single FF (or not — the honest result).

---

## Idea 03 — Batched ColabFold / AlphaFold Proteome Scan
*Mining an uncharacterized protein family for folds and functions.*

- **Question:** Predict high-confidence structures for every protein in an understudied family
  or pathogen proteome and infer function via fold-matching.
- **Method:** One ColabFold (AlphaFold2) prediction per protein as a GPU array task; score by
  pLDDT; annotate via Foldseek structural search against PDB100.
- **Cluster role:** Hundreds of proteins × 15–45 min GPU each. AlphaFold is embarrassingly
  parallel *at the protein level* — precisely what `qsub -t` does. ~1–2 days on 8 GPUs.
- **Originality:** Applies cutting-edge structure prediction to a *family*, not a single
  target, converting "structure prediction" from a toy into a functional-annotation pipeline.
- **Deliverable:** Predicted-structure gallery, fold-match/function table, family clustering
  dendrogram.

---

## Idea 04 — Virtual Drug-Repurposing Screen
*Ranked repurposing shortlist from docking, with honesty about what docking can't tell you.*

- **Question:** Find approved drugs that dock well against a disease target, as a cheap
  computational repurposing screen.
- **Method:** Dock the ~3,000-compound DrugBank approved set against a target pocket with
  AutoDock Vina (CPU) or GNINA (GPU); calibrate the score cutoff on positive controls; cluster
  hits by scaffold; re-rank the top ~50 at higher exhaustiveness.
- **Cluster role:** 1 task per ligand; ~500–1,500 CPU-hours or a few GPU-hours compressed into
  ~1 day on the cluster.
- **Originality:** The scientific content is the *honest framing* — a ranked, chemotype-diverse
  shortlist with a calibrated cutoff and an explicit "docking scores are not binding free
  energies" caveat — which is exactly the rigor most repurposing fan-projects lack.
- **Deliverable:** Ranked repurposing table + score histogram with positive-control threshold +
  3D pose figures.

---

## Idea 05 — Deep-Mutational-Scan Validation
*A rigorous benchmark of empirical predictors against a published, free experimental dataset.*

- **Question:** How well do FoldX/Rosetta reproduce a published experimental deep mutational
  scan (~1,000 mutations), and where do they fail?
- **Method:** Predict every mutation in a DMS dataset (e.g., TEM-1, GB1, GFP, SUMO) with FoldX
  and Rosetta `cartesian_ddg`; compute Pearson, Spearman, and ROC against experiment; dissect
  error by burial, substitution class, and secondary structure.
- **Cluster role:** ~1,000 × 2 predictor runs ≈ 80–160 CPU-hours → hours on an array.
- **Originality:** It is a *methods-validation* project — the same shape as a published
  benchmark — and the error-anatomy (where predictors break) is more informative than the
  headline correlation.
- **Deliverable:** Predicted-vs-experimental correlation/ROC figure + error-anatomy table.

---

## Idea 06 — Binding Free-Energy Consensus (MM/PBSA + TI/MBAR)
*Two method families, one set of ligands: which ranks binders better?*

- **Question:** Compare rigorous alchemical TI/MBAR binding free energies against the cheap
  end-state MM/PBSA method on the same protein–ligand complexes, and against experiment.
- **Method:** Build 5–8 well-characterized complexes; run both methods (with replicates and
  force fields); compare rank order and absolute accuracy.
- **Cluster role:** TI/MBAR is ~200 GPU-window-hours per ligand set; MM/PBSA is a CPU array.
  The two-method × multi-replicate design is only affordable as an array.
- **Originality:** The comparison *between* method families (not just predictor vs
  experiment) is the contribution — cheap-method agreement with rigorous-method is a decision
  tool for anyone running virtual screening.
- **Deliverable:** Predicted vs experimental ΔG (both methods, error bars) + rank-agreement
  verdict.

---

## Idea 07 — Peptide-Library Interface Scan
*In-silico screening of a peptide library for competitive binders at a protein interface.*

- **Question:** Which short peptides bind a therapeutically relevant protein–protein
  interface (e.g., SOD1 dimer interface, ACE2–spike, a PDZ domain)?
- **Method:** Two-stage: dock a systematic peptide library (2,000–8,000) with Vina/HADDOCK,
  then validate the top ~20 with short MD or MM/PBSA.
- **Cluster role:** Thousands of small docking jobs (the wide stage) plus a handful of
  expensive validation runs (the narrow stage) — a two-stage design that only works because
  the wide stage ran in parallel.
- **Originality:** The systematic enumeration (e.g., all single-point variants of a known
  binder) yields a per-position preference matrix — a mini interface-level saturation scan.
- **Deliverable:** Ranked library + validation table + per-position design rule.

---

## Idea 08 — Coarse-Grained Self-Assembly Sweep
*A phase diagram of mesoscale self-assembly from MARTINI MD.*

- **Question:** How does lipid/peptide/surfactant self-assembly depend on composition,
  temperature, and seed?
- **Method:** Coarse-grained GROMACS + MARTINI simulations across a composition × temperature
  grid, multiple seeds each; classify final morphologies (micelle / vesicle / bilayer /
  fiber).
- **Cluster role:** 1 array task per (condition, seed); CG MD is 100× cheaper than atomistic
  but still hours each, and 5 seeds per condition is the difference between noise and a phase
  diagram.
- **Originality:** A genuine 2-D phase diagram as the headline result — the grid is the
  science, and the grid requires a cluster. Ties to the aggregation biology of the parent
  project.
- **Deliverable:** Phase diagram (composition × temperature → morphology) + kinetics +
  morphology gallery.

---

## Idea 09 — Kinetic Stability from Unfolding Replicates
*Time-to-first-unfold as a disease-relevant metric beyond ΔΔG.*

- **Question:** Which variants unfold fastest, not just which are least stable? Kinetic
  stability may be the more disease-relevant quantity for aggregation-prone proteins.
- **Method:** Replicate accelerated-unfolding MD (elevated temperature, order-parameter
  thresholding on Q-fraction); Kaplan–Meier-style survival curves; mean first-passage times;
  correlate with thermodynamic ΔΔG.
- **Cluster role:** 10 variants × 10+ replicates of 200–500 ns GPU runs = the only way to get
  statistically meaningful rare-event survival curves.
- **Originality:** Most projects measure *depth* of the well (ΔΔG); this measures *escape
  rate*, and the kinetic-vs-thermodynamic divergence cases are the scientific gold.
- **Deliverable:** Survival curves per variant + kinetic-vs-thermodynamic comparison table.

---

## Idea 10 — Dimerization Free-Energy Panel
*Extending FEP to the dimer interface — the hardest physics in the series.*

- **Question:** Do ALS mutations destabilize the SOD1 dimer interface, and can alchemical FEP
  quantify ΔΔG of dimerization?
- **Method:** Alchemical mutation in both dimer and monomer systems, with a correctly handled
  dissociation reference and symmetry factor; ~10 variants × 2 systems × 2 legs × 18 windows ×
  3 reps.
- **Cluster role:** ~2,000+ GPU windows; the entropically hard dissociation leg needs
  extensive sampling. The most expensive project here.
- **Originality:** Directly extends the parent repo into new physics (interface stability) and
  connects to Lindberg's class-1/class-2 framework — producing a 2-D monomer × dimer ΔΔG map
  that classifies variants mechanistically.
- **Deliverable:** ΔΔG_dimerization table vs experiment + interface-vs-buried scatter +
  variant classification map.

---

## Idea 11 — Thermostability Redesign
*Consensus FEP as a protein-engineering tool — designing stabilizing mutations, not finding
destabilizing ones.*

- **Question:** Can consensus-FF alchemical FEP rank candidate mutations by *stabilizing*
  potential, validated against a known literature-stabilizing mutation?
- **Method:** FoldX pre-screen → FEP the top ~20 candidates (2 FFs × 2 legs × 18 windows × 3
  reps) → rank; validate the protocol on one published stabilizing mutation.
- **Cluster role:** ~4,000 GPU windows for 20 candidates; the consensus requirement makes it
  2–3× the single-FF cost but is essential for credibility.
- **Originality:** Flips the disease-variant question into an engineering question; a ranked,
  error-barred stabilizing-mutation shortlist is exactly what enzyme engineers need.
- **Deliverable:** Ranked stabilization shortlist + FEP-vs-FoldX agreement + protocol
  validation on the literature mutation.

---

## Idea 12 — Cryptic-Pocket Detection from an MD Ensemble
*Finding transient druggable pockets the crystal structure hides.*

- **Question:** Do transient "cryptic" binding pockets open on a protein that looks
  undruggable in its crystal structure?
- **Method:** Many independent MD replicas (32 × 1 µs); fpocket on every frame; cluster into
  persistent pocket identities; rank by open-frequency × residence × volume; stretch: dock a
  fragment library into the top pocket.
- **Cluster role:** ~30 GPU-days of atomistic MD; the entire premise is *many* long replicas
  because rare events need ensemble statistics.
- **Originality:** Applies the modern "undruggable target" narrative to a real protein (SOD1's
  aggregation-prone loops are a natural target) and requires a positive control (a known
  cryptic-pocket protein) to validate the pipeline.
- **Deliverable:** Pocket-volatility map + ranked cryptic-pockets table + fragment-docking
  evidence.

---

## Idea 13 — Host–Pathogen Protein-Protein Docking Scan
*A systematic interaction-candidate screen between a pathogen's proteome and host proteins.*

- **Question:** Which pathogen proteins plausibly interact with human proteins, and what
  pathways do they target?
- **Method:** ColabFold all pathogen proteins (Idea 03), then HADDOCK/ClusPro every
  (pathogen × host) pair; calibrate on a known interaction (e.g., spike–ACE2); cluster top hits
  by host function.
- **Cluster role:** ~500 pairs × 1–3 h CPU each = ~1,000 CPU-hours → days on an array.
- **Originality:** Converts hypothesis-generation into a systematic screen with a calibrated
  cutoff and functional enrichment — hypothesis *generation with statistics*, not a list of
  guesses.
- **Deliverable:** Ranked interaction table + positive-control calibration + host-pathway
  enrichment.

---

## Idea 14 — ML Surrogate Trained on Cluster-Generated Stability Data
*Can a model learn the output of expensive physics and predict instantly?*

- **Question:** Train an ML surrogate on thousands of cluster-generated ΔΔG values and
  benchmark it against the physics that generated them, plus a physics+ML hybrid.
- **Method:** Generate labels via saturation scan (cheap) and/or FEP (gold); engineer sequence
  and structural features; grid-search boosting and/or a GNN; split by protein (no leakage).
- **Cluster role:** The *dataset is generated by arrays*; the hyperparameter sweep is another
  array. The whole "generate → train → evaluate" loop is cluster-scale.
- **Originality:** Most student ML projects use someone else's dataset; here the student
  *creates* the dataset with a cluster, and the physics-vs-ML-vs-hybrid comparison is the
  result — a learning curve showing the value of the compute they control.
- **Deliverable:** Surrogate-vs-physics-vs-hybrid benchmark + learning curve + fast scorer.

---

## Idea 15 — Soft-Core Protocol Benchmark
*Which alchemical soft-core scheme is more accurate and robust? A methods decision, with
numbers.*

- **Question:** Beutler vs Gapsys soft-core: which gives more accurate ΔΔG and converges more
  robustly for protein mutations?
- **Method:** Identical FEP windows differing only in the soft-core `.mdp` parameters; measure
  accuracy vs experiment and robustness (hysteresis, overlap, end-state stability).
- **Cluster role:** 12 mutations × 2 schemes × 2 legs × 18 windows × 3 reps ≈ 2,600 GPU
  windows — a controlled experiment that is only meaningful because every arm is run at equal
  cost and scale.
- **Originality:** A protocol-recommendation paper for the parent repo (and the field); the
  least flashy and most publishable. Directly answers a question the parent config flags as
  open (pipeline.yaml:116-127).
- **Deliverable:** Beutler-vs-Gapsys accuracy and robustness comparison + a written config
  change proposal.

---

# Series II — The Biophysics × Data-Science Series (16–30)

## Idea 16 — Conformal Prediction for ΔΔG
*Error bars with a guaranteed coverage promise — the honest-uncertainty project.*

- **Question:** Can distribution-free conformal prediction wrap FoldX (or FEP) ΔΔG predictions
  into calibrated, *guaranteed-coverage* intervals — and do those agree with MBAR/replicate
  errors?
- **Method:** Calibrate a conformal interval on a large FoldX scan (split by position); measure
  empirical coverage overall and in mutation-class slices; compare widths to MBAR errors on the
  FEP subset.
- **Cluster role:** The calibration set is a ~1,000–2,000 mutation array; conformality without
  scale is meaningless.
- **Originality:** Every FEP pipeline reports an error bar, but almost nobody *audits* whether
  it covers the truth. Conformal prediction is a modern ML guarantee few students know, and the
  "where do conformal and MBAR disagree" result is genuinely useful.
- **Deliverable:** Coverage-vs-nominal calibration plot + interval-width comparison table.

---

## Idea 17 — Active Learning for the Validation Gate
*Which controls should get FEP'd first? Spend GPU-hours where they inform the decision most.*

- **Question:** Given the go/no-go gate, which next variant maximizes gate information per
  GPU-hour?
- **Method:** Fit a cheap surrogate over all controls; query variants by uncertainty sampling /
  expected gate-metric change; run real FEP on the chosen few; update; repeat; compare against
  random and fixed-set baselines.
- **Cluster role:** Each round is a batch of FEP arrays; evaluating multiple query strategies
  multiplies the compute — only a cluster affords the loop.
- **Originality:** Applies active-learning / Bayesian experimental design to a *decision a real
  pipeline has to make*, producing a power-style "how many controls does the gate really need"
  conclusion.
- **Deliverable:** Gate-metric-vs-GPU-hours learning curves per strategy + query map.

---

## Idea 18 — Physics-Informed ML with Thermodynamic Cycle Consistency
*A stability predictor that must respect the identity ΔΔG(A→B) = −ΔΔG(B→A).*

- **Question:** Does enforcing thermodynamic cycle consistency (soft penalty or hard
  node-potential parameterization) improve a stability predictor, especially on extrapolation?
- **Method:** Build a mutation-graph predictor; compare unconstrained vs consistency-penalized
  vs hard-consistent versions on held-out positions; quantify cycle-residual imbalance.
- **Cluster role:** Labels come from saturation-scan arrays; cycle/triple enumeration and the
  training sweep are parallel.
- **Originality:** The "your predictor violates thermodynamics by this much" diagnostic is a
  memorable, defensible contribution — physics-informed ML is a frontier topic in 2026 ML.
- **Deliverable:** Cycle-imbalance diagnostic + constrained-vs-baseline accuracy comparison.

---

## Idea 19 — Markov State Models of the SOD1 Folding Landscape
*Turn many short trajectories into a kinetic network — free energies, pathways, and rates.*

- **Question:** Do disease variants change the *kinetic network* of SOD1 — populations,
  barriers, transition paths — in ways single-trajectory analysis misses?
- **Method:** Build MSMs (PyEMMA) from 60+ short trajectories per variant: featurize → tICA →
  cluster → MSM → validate (implied timescales, Chapman–Kolmogorov) → extract stationary
  distribution, landscapes, mean first-passage times.
- **Cluster role:** MSMs *require* many independent trajectories; the array of 240 GPU runs is
  the data-generator. One long run is strictly worse than many short ones.
- **Originality:** Most student MD is "one trajectory, one story." MSMs are ensemble kinetics —
  the statistically honest way to describe a landscape, and a genuinely advanced method.
- **Deliverable:** Per-variant free-energy landscapes + kinetic-network diagrams + rate
  comparison table.

---

## Idea 20 — Learned Collective Variables with VAMPnets
*Let a neural network discover the slow coordinates of protein motion.*

- **Question:** Can a time-lagged autoencoder (VAMPnet) learn collective variables that beat
  hand-crafted ones (RMSD, Q, tICA) at separating metastable states?
- **Method:** Train a VAMPnet on the same trajectory ensemble as Idea 19; score with VAMP-2 on
  held-out trajectories; compare landscapes and add saliency interpretation.
- **Cluster role:** Shares Idea 19's trajectory array (one data-generating run, two analyses);
  training sweep is another array.
- **Originality:** Representation learning for molecular dynamics is cutting-edge; a student
  who can explain *why* a learned CV beats Q is already ahead of most undergraduates.
- **Deliverable:** VAMP-2 comparison + landscape projections + interpretability (saliency)
  pass.

---

## Idea 21 — Epistasis Detection: When Is a Double Mutant Not the Sum of Its Parts?
*The assumption that single-mutant scans build on, tested systematically.*

- **Question:** Can we detect non-additive coupling (epistasis) between mutation pairs, and
  does structural distance predict it?
- **Method:** Fit ΔΔG ≈ Σ singles + Σ pair-terms over single AND double mutants (FoldX),
  generating the ~N² double-mutant matrix by array; correlate pair-epistasis with structural
  features; validate against experimental double-mutant DMS (GB1).
- **Cluster role:** ~5,000–20,000 double-mutant tasks — an N² blow-up only a cluster absorbs.
- **Originality:** Every triage tool assumes additivity; this project measures *where that
  assumption breaks* — a foundational, under-appreciated question.
- **Deliverable:** Epistasis heatmap + strongest-pairs table + additivity-failure rule.

---

## Idea 22 — Protein-Language-Model Embeddings as Stability Features
*The "pLM vs physics" debate, tested on your data.*

- **Question:** Do ESM-2 zero-shot scores predict stability as well as structure-based physics,
  and is the combination better than either?
- **Method:** Score a full saturation scan with ESM-2 log-likelihoods; compare physics-only,
  pLM-only, and hybrid predictors on held-out mutations; dissect where each wins.
- **Cluster role:** ESM-2 scoring over a protein family is thousands of GPU tasks; model sweeps
  are another array.
- **Originality:** The hybrid test — are the signals *redundant or complementary*? — is a
  genuinely argued literature question, and the per-class error anatomy gives a real answer.
- **Deliverable:** Single-source-vs-hybrid benchmark + error-anatomy table.

---

## Idea 23 — The Predictability Ceiling
*Even a perfect predictor can't beat experimental noise — how much is knowable at all?*

- **Question:** What maximum correlation/MAE is achievable given the noise floor of the
  experimental ΔΔG measurements themselves?
- **Method:** Estimate σ_exp from repeated measurements of the same mutations (ProTherm,
  DMS replicates, the parent repo's multi-source controls); simulate the noise-corrected
  ceiling; score real predictors against it as a "fraction achieved."
- **Cluster role:** Thousands of bootstrap/Monte-Carlo resamples — an array of cheap,
  independent tasks.
- **Originality:** A meta-science / information-theory contribution: it reframes how the parent
  repo's gate (and the whole field) should read predictor performance. Few students even know
  the concept of a noise floor exists.
- **Deliverable:** The ceiling figure (max achievable R vs σ_exp, predictors plotted against
  it) + a corrected performance table.

---

## Idea 24 — Transfer Learning Across Proteins
*Can one protein's cluster-generated physics teach a model about another?*

- **Question:** For a protein with few experimental points, is it better to train only on those,
  transfer from another protein's labeled scan, or fine-tune?
- **Method:** Generate FoldX scans for 2–4 source proteins; evaluate three regimes
  (no-transfer / direct-transfer / fine-tune) across many splits and target proteins; relate
  transfer benefit to source-target structural similarity.
- **Cluster role:** Source scans are arrays; the train/test grid is hundreds of independent
  runs.
- **Originality:** The *transfer curve* (error vs target-data-points per regime) plus the
  structural-similarity hypothesis is a real, decision-relevant result for anyone doing
  small-data protein prediction.
- **Deliverable:** Transfer curves + benefit-vs-similarity scatter + a practical labeling
  recommendation.

---

## Idea 25 — Are FEP Error Bars Honest? A Bootstrap Audit of MBAR
*The statistics-of-simulation project: do the reported uncertainties actually cover the
truth?*

- **Question:** Compare MBAR's analytic error to block-bootstrap, replicate-bootstrap, and
  chain-bootstrap errors on real FEP data — and find where they diverge.
- **Method:** Run a real FEP panel with ≥5 replicates for a subset; refit MBAR thousands of
  times under each resampling scheme; audit per-variant error estimates.
- **Cluster role:** Extra replicates are real GPU work; the ~10⁴ bootstrap refits are an
  embarrassingly-parallel array.
- **Originality:** A calibration report ("your error bars are honest/optimistic by factor X")
  that directly answers the parent repo's own rule-5 mandate — and a candidate fix (raise the
  replicate floor).
- **Deliverable:** Error-estimator audit table + bias-vs-variance decomposition +
  replicate-floor recommendation.

---

## Idea 26 — Predicting FEP Window Failure Before You Pay for It
*Reliability ML: skip the windows that are going to die.*

- **Question:** Can a classifier predict which (variant, leg, window, replicate) will crash or
  fail to converge — from cheap pre-run features including a 10 ps probe?
- **Method:** Run a real panel while logging failures; extract cheap features (mutation
  context, λ, soft-core, probe diagnostics); train a failure classifier; test mitigations
  (rescue rate).
- **Cluster role:** Label-generating array + a 100×-cheaper probe array + training grid; the
  "prediction + mitigation" loop only pays off at cluster scale.
- **Originality:** Turns a nuisance (failed windows) into a modeling problem with a concrete
  savings payoff — an applied-reliability angle almost no student project touches.
- **Deliverable:** Failure-prediction ROC/precision + feature-importance + rescue-rate result.

---

## Idea 27 — Co-Evolution Coupling vs FEP
*Do the free sequence signal and the expensive physics agree — and when they disagree, why?*

- **Question:** Do DCA/EVcouplings evolutionary couplings predict destabilization, and how does
  that compare with FEP/FoldX on the same variants?
- **Method:** Deep MSA → coupling matrix; run cheap (FoldX, DCA) and expensive (FEP) signals on
  a variant panel; correlate all three against experiment; structurally dissect the
  disagreements.
- **Cluster role:** The FEP panel (2,200+ windows) is the cost; DCA is cheap; the array is the
  only way to afford the panel at a meaningful size.
- **Originality:** Multi-signal fusion with an explicit *tension analysis* — "DCA says coupled,
  FEP says neutral: functional constraint or missing physics?" — plus a practical conclusion
  (can DCA cheaply pre-screen VUS before FEP?).
- **Deliverable:** Three-signal comparison + disagreement case studies + pre-screen
  recommendation.

---

## Idea 28 — Gaussian-Process Surrogates for the λ-Ladder
*Adaptive λ-placement: spend windows where the free energy is nonlinear.*

- **Question:** Can a GP model of the running windows tell you adaptively where to add λ
  windows, saving GPU time at equal (or better) accuracy than the fixed 18-window ladder?
- **Method:** GP over dG(λ) (and u_kn overlap); acquisition function picks the next λ; repeat
  over variants; compare against fixed ladder at equal window count.
- **Cluster role:** Each adaptive episode is a sequence of window arrays; the repeat study
  across variants multiplies the compute.
- **Originality:** Active learning applied *inside* the FEP estimator — a true methodology
  contribution, and the mock harness in the parent repo makes it developable cheaply.
- **Deliverable:** Fixed-vs-adaptive error-vs-windows comparison + learned sampling pattern +
  config recommendation.

---

## Idea 29 — Information-Theoretic Frame Selection
*Which frames actually matter for MBAR — and how little can you store?*

- **Question:** Which frame-selection strategy minimizes MBAR variance for a fixed storage
  budget?
- **Method:** Compare decorrelation, thinning, variance-minimizing selection, and u_kn-space
  clustering on real long windows at several budgets; measure ΔΔG error vs frames kept.
- **Cluster role:** Long real windows plus ~10⁴ parallel MBAR refits; directly informs the
  repo's `trajectory_retention` policy.
- **Originality:** Most pipelines thin frames arbitrarily; the information-theoretic answer
  ("keep X% selected by Y") is an actionable, disk-saving result with a clean error trade-off
  figure.
- **Deliverable:** Error-vs-frames-kept curves + a concrete retention recipe.

---

## Idea 30 — Anomaly Detection on MD Trajectories
*Let unsupervised ML find the rare events standard order parameters wash out.*

- **Question:** Can autoencoder/one-class anomaly detection flag disease-relevant rare
  conformational events that RMSF/Q miss — and do variants change the anomaly budget?
- **Method:** Train an anomaly model on WT frames only; score all frames of WT + variant
  trajectories; compare anomaly rates, durations, and per-residue localization; validate
  against a physical order parameter.
- **Cluster role:** 120–240 GPU trajectories (the rare events need ensemble statistics) plus a
  parallel anomaly-scoring pass per trajectory.
- **Originality:** Framing MD analysis as anomaly detection (what does "abnormal" look like
  when "normal" is learned, not defined) is a modern, memorable angle, and the WT-only training
  discipline is a subtle design point that impresses.
- **Deliverable:** Anomaly-rate comparison (WT vs variants) + per-residue anomaly map + top-
  event case studies.

---

# How the two series compose

The strongest project is often a *combination* rather than a single idea. The README lists the
highest-value pairings; the most natural thesis-shaped ones:

| Pairing | Why it works |
|---------|--------------|
| **25 + 02** | Bootstrap-audit the consensus-FF error bars — methods rigor on your own data. |
| **16 + 01** | Conformal intervals around a saturation scan — guaranteed-coverage stability maps. |
| **19 + 20** | One trajectory array feeds both an MSM and a VAMPnet — kinetic networks plus learned coordinates. |
| **28 + 02** | Adaptive λ-placement evaluated *inside* the consensus-FF benchmark — two methods questions in one run. |
| **26 + 10** | Predict dimer-FEP failures before paying for the most expensive array in the series. |
| **22 + 05** | pLM embeddings vs physics on a deep mutational scan — the modern debate, tested. |
| **17 + 02** | Active learning picks which controls get consensus-FF first — optimal gate design. |

---

*Prepared as a planning digest. For full research questions, fan-out geometry, cluster
budgets, milestones, deliverables, and pitfalls, see the individual `idea-NN_*.md` files and
the `README.md` matrix/ladder.*
