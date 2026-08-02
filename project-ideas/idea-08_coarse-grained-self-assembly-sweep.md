# Idea 08 — Coarse-Grained Self-Assembly Sweep

## Research question

**How does self-assembly depend on composition, concentration, and temperature?** Using
coarse-grained (CG) molecular dynamics — where groups of atoms become single beads — you can
simulate thousands of molecules for microseconds and watch them assemble into micelles,
vesicles, bilayers, fibers, or clusters. A parameter sweep (composition × temperature ×
concentration, multiple seeds) maps the **phase/assembly behavior** of your system.

Framing options:
- **Lipid composition → vesicle/bicelle morphology** (how does tail length / headgroup ratio
  change what forms?).
- **Peptide amphiphiles / surfactant mixtures → micelle vs fiber** (relevant to the parent
  project's amyloid angle!).
- **Model protein CG beads → aggregate vs remain dispersed** as a function of a "hydrophobicity"
  parameter — a toy model of protein aggregation (ties directly to SOD1 aggregation thinking).

## Why this needs an SGE cluster

```
tasks = compositions × temperatures × seeds  (tens to hundreds)
```

- **1 array task = 1 (system, seed) run.** CG MD is ~100× cheaper than atomistic, but each
  run is still hours of GPU (or threaded CPU) — and you want **many seeds** because assembly
  is stochastic: 5 seeds per condition is the difference between a noise artifact and a phase
  diagram.
- The output is a **phase diagram** — a genuinely 2-D scientific result that only exists if
  you ran the whole grid, and the grid is the cluster's job.
- CG GROMACS runs well on `-pe omp` CPU nodes too, so this can use otherwise-idle cores
  alongside the GPU jobs.

## Data & tools

- **Tools:** GROMACS + **MARTINI** CG force field (free; `martinize.py` builds topologies;
  MARTINI 3 is current). Vesicle/bilayer building via `insane.py` or CHARMM-GUI MARTINI.
- **Systems:** DPPC/DOPC/cholesterol lipid mixtures; single-tail surfactant; CG peptides
  (MARTINI protein beads); or fully custom bead "designer amphiphiles".
- **Analysis:** cluster detection (`gmx cluster`, or MDAnalysis DBSCAN), shape metrics
  (radius of gyration, ellipticity, radial density), size distributions over time.

## Skill prerequisites

- Basic MD literacy (what a trajectory and a topology are).
- Python (MDAnalysis/matplotlib) for the phase-diagram analysis.
- Willing to learn MARTINI conventions (bead types, elastic network).

## Cluster budget

| Parameter | Value |
|---|---|
| Conditions | 16 (4 compositions × 4 temps) |
| Seeds per condition | 5 |
| Runs | 80 × ~4–8 h GPU (or 8–16 h CPU-threaded) |
| Wall-clock on 8 GPUs | **~2–4 days** |

## Milestones

1. Build one system (e.g., a lipid bilayer patch) with MARTINI; equilibrate; verify it
   behaves (area per lipid in range).
2. Define the sweep grid (composition × temperature); for each, pick a starting box and
   N molecules.
3. Run 3 seeds of one condition; check reproducibility (do all 3 assemble to the same
   morphology?).
4. Array over (condition, seed); each writes trajectory + a per-frame structure metric
   (cluster count/size) so the "estimate" is small.
5. Aggregate into a **phase diagram**: color each (composition, temperature) cell by the
   dominant structure (micelle / vesicle / bilayer / fiber / dispersed).
6. Characterize transitions: does a composition change flip the morphology? A temperature
   change? How sharp?

## Deliverables

- **Phase diagram** (composition × temperature → structure) — the headline figure.
- Morphology gallery (representative snapshots per phase).
- A time-course showing assembly kinetics (size vs time per condition).
- A short interpretation: what drives the phase boundary?

## Pitfalls

- **Seed dependence is real.** Never draw conclusions from one trajectory; 5 seeds minimum,
  and plot *distributions*, not single runs.
- **MARTINI is not atomistic.** It cannot capture hydrogen bonds or specific electrostatics;
  it is for *mesoscale morphology*. Say so explicitly.
- **Finite-size effects.** A box of 200 molecules can't form a giant vesicle. Report the box
  size and how morphology depends on it.
- **Equilibration is slow in CG too.** Assembly often takes µs of CG time; a run that ends at
  a metastable micelle is a real result, not a failure — but you must know it's metastable
  (check the cluster size plateaus).
- **Checkpoint.** Long CG runs WILL be preempted at the wall-time limit; write checkpoints or
  split into extendable legs (`-cpi`).
