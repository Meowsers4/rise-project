# F64A selection-policy sensitivity results — 2026-09-16

## Evidence and integrity

The user completed the authorized CPU analysis on SCC and supplied the full
[report](f64a_extension_review/f64a_selection_sensitivity_v1.json). Its SHA-256 is
`10490de20f4f075b7db66f6c18a15611285e482db29525ca1a86ad7077aa3866`.
Analysis commit: `901bfb313dbf3d679c10a4010960692b0a02e0e0`;
simulation commit: `80097dcb2ff90ebfdedfadb6e4d9a866a3d9b55a`;
simulation protocol: `bf6841ccb3b9de79`.

All 18 fits, 360 window-selection records, 342 adjacent-pair records and 171 paired
local-change records are present. Local review verified the seven recorded source/
config/reference hashes against the repository, all 60 NPZ input hashes against
the frozen block report, paired mean/SEM arithmetic, and all six adaptive reference
checks including exact window diagnostics. The report records clean analysis paths
and successful adaptive reproduction. This is review of a user-supplied SCC report,
not independent local MBAR reproduction; the extended matrices remain on SCC.

The [design](f64a_selection_sensitivity_design.md) was fixed before this CPU run,
after viewing the adaptive block results. It is exploratory, not blind validation.

## Result

Compare retained raw columns `[6001:7501]` and `[7501:9001]`, without duplicating
the boundary. All entries below are folded-leg mutation free-energy movements
(later minus earlier), in kcal/mol; they are not folding ΔΔGs.

| Repeat | Adaptive trimming | No additional trimming | Fixed 25% trimming |
|---|---:|---:|---:|
| r0 | −0.352654 | +0.158341 | −0.109051 |
| r1 | −0.565827 | −0.836415 | −0.694710 |
| r2 | −0.729934 | −0.776855 | −0.725439 |
| Mean | −0.549472 | −0.484976 | −0.509733 |
| Paired repeat SEM | 0.109218 | 0.322118 | 0.200537 |

The mean movement spans only 0.064496 kcal/mol across the three policies. Downward
movement is therefore not specific to adaptive cutoff selection. However, **not
every repeat moves downward under every policy**: r0 changes sign without trimming.
r1 and r2 move downward under all three policies. Selection sensitivity is reduced
as an explanation of the mean, not eliminated as a concern for individual repeats.

Paired SEMs describe variability across three repeat movements; they are not 95%
confidence intervals or independent-sample significance tests. The report also
contains each fit's conditional MBAR uncertainty. Policies share the same raw data,
and disjoint time blocks can remain correlated through slow modes.

## Sample selection and local diagnostics

Fixed policies do not necessarily retain more samples after correlation thinning.
Their smallest per-window selected counts range from 13 to 26 across fits, and
maximum estimated statistical inefficiencies reach 124.1 without trimming and
81.4 with a fixed 25% trim. Window 18 repeatedly has low selected counts, but other
windows (including 9, 12, 13, 14 and 17) also matter. Minimum adjacent overlap spans
0.0126–0.0644 across all policies/fits. No solver notes were recorded; that does not
certify physical equilibration or reliable correlation estimates on changing data.

Adjacent Zwanzig changes do not identify a unique, policy-independent target:
the largest pooled mean absolute bidirectional changes are at pair 18–19 under
adaptive selection, and at 14–15 under both fixed policies. Individual repeats
also differ. Forward and reverse changes can oppose each other. These are noisy
energy-space diagnostics, not contributions that sum to the MBAR endpoint drift,
proof of a causal bottleneck, or evidence for a particular protein conformational mode.

## Interpretation and proposed next step

The mean late movement persists across the tested cutoff policies, and two repeats
show consistent downward movement. Duration-stable estimates at 9 ns remain
unestablished. Small conditional uncertainties cannot account for configurations
never visited, as emphasized in the
[alchemical best-practices guide](https://pmc.ncbi.nlm.nih.gov/articles/PMC8388617/).
These data do not isolate physical slow sampling from every analysis sensitivity,
nor do they diagnose force-field or experimental-reference error.

The recommended next test is a bounded **9→15 ns retained-sampling continuation of
the full F64A folded ladder**, all 20 windows and three repeats, rather than a
selective extension based on local hysteresis. The
[draft design](f64a_continuation_9_to_15_design.md) is for user consideration;
documentation/design authorization is not GPU submission authorization.
The historical failed gate, F64A exclusion and apo-2SH state are unchanged.

Per the user's current instruction, Peter Vlasov's review is deferred until the
user reintroduces it. It is not a prerequisite or pending action in this workflow;
no reviewer has been contacted.
