# Reproducible forensic figure package

Status: complete, CPU-only, 2026-09-14. This package builds Figures 1--3 and 5 in the
methods-and-limitations outline directly from archived JSON and topology files. Figure 4 is
the separately executable reference-sensitivity analysis. None of these figures is a new gate
evaluation, and the gate-of-record verdict remains failed.

## Outputs

| figure | output | source evidence |
|---|---|---|
| 1 | [`gate_of_record.png`](figures/gate_of_record.png) | frozen attempt-1 `ddg.json` records; F64A exclusion re-applied from the registered config |
| 2 | [`convergence_diagnostic.png`](figures/convergence_diagnostic.png) | 24 folded-leg records in independently re-derived convergence JSON |
| 3 | [`a4v_sampling_sensitivity.png`](figures/a4v_sampling_sensitivity.png) | archived A4V 3 ns result versus retained A4V 9 ns result |
| 4 | [`reference_sensitivity.png`](figures/reference_sensitivity.png) | audited reference scenarios; see [`reference_sensitivity_analysis.md`](reference_sensitivity_analysis.md) |
| 5 | [`g93a_disulfide_diagnostic.png`](figures/g93a_disulfide_diagnostic.png) | frozen 2SH and SS convergence JSON plus all six folded `hybrid.top` files |

Figure 2 computes each point as
`abs(DG_i - mean(DG_j for the two independently solvated siblings))`. It reproduces Pearson
`+0.071866` and Spearman `+0.080000` over 24 folded records. Figure 5 parses the atom and bond
tables rather than trusting residue labels: all three 2SH topologies lack the target SG--SG
bond and retain both target HG atoms; all three SS topologies contain the bond and lack both
HG atoms. The folded mean shift is `-0.016918 kcal/mol`.

## Rebuild

From the repository root, with the two local archives mounted at their recorded locations:

```bash
python -m src.analysis.methods_figures \
  --attempt1-archive /Users/bodebosell/sod1fep_archive_2026-09-11 \
  --reanalysis-convergence /Users/bodebosell/sod1fep_archive_2026-09-11/reanalysis_2026-09-12/convergence \
  --ss-archive /Users/bodebosell/sod1fep_archive_2026-09-14/G93A_SS_diagnostic \
  --output-dir docs/figures
```

The builder stops rather than plotting if the frozen archive no longer reproduces the failed
seven-point gate, the folded diagnostic does not contain exactly 24 records, either G93A arm
lacks three folded replicas, or the topology evidence does not match the declared redox arm.

## Archive integrity

[`forensic_archive_manifest.json`](../data/forensic_archive_manifest.json) records SHA-256
checksums for every raw NPZ retained in the main archive and SS diagnostic, plus the direct JSON,
topology, configuration, and builder inputs. The main archive has 840 re-derivable attempt-1
windows for the seven non-A4V variants and 120 A4V attempt-2 windows; A4V's attempt-1 windows
were overwritten and are not misrepresented as recoverable. The SS archive adds 120 windows.
Paths are stored relative to three named roots, so the manifest does not depend on one machine's
absolute path. The aggregate checksum covers each relative path, file size, and individual
checksum.
