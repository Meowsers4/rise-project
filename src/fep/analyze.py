"""Stage 3 (analysis) -- windows -> per-variant ΔΔG via MBAR, with uncertainty.

Aggregates the per-window reduced potentials (from :mod:`src.fep.window`) over both
legs and all replicates, runs MBAR per (leg, replicate), and forms

    ΔΔG_folding = ΔG(folded leg) - ΔG(unfolded leg)      (positive = destabilizing)

averaged over replicates. A ΔΔG is reported only WITH an uncertainty and a
convergence check (CLAUDE.md #5): the forward/reverse (Zwanzig) hysteresis of each
leg is the available per-variant convergence proxy, and ``converged`` is False if it
exceeds ``fep.convergence.max_cycle_closure_kcal``.

Emits ``results/fep/<variant>/ddg.json`` in the schema the validation gate consumes
(see :data:`src.analysis.validate.FEP_DDG_SCHEMA`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymbar import MBAR

from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]
_KB_KCAL = 0.0019872041  # Boltzmann constant, kcal/mol/K


def _logmeanexp(a: np.ndarray) -> float:
    """log(mean(exp(a))) with the max-shift trick (avoids overflow)."""
    m = np.max(a)
    return float(m + np.log(np.mean(np.exp(a - m))))


def load_leg_replicate(fep_dir: Path, variant: str, leg: str, rep: int,
                       n_states: int) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Aggregate a leg's windows for one replicate into MBAR inputs.

    Returns (u_kn, N_k, per_window_u) where u_kn is (n_states, total_samples) and
    per_window_u[k] is window k's (n_states, n_k) block (used for hysteresis).
    """
    per_window = []
    for w in range(n_states):
        npz = np.load(fep_dir / variant / leg / f"w{w}_r{rep}.npz")
        u = npz["u_kn_window"]
        assert u.shape[0] == n_states, f"{variant}/{leg}/w{w}_r{rep}: state mismatch"
        per_window.append(u)
    N_k = np.array([u.shape[1] for u in per_window])
    u_kn = np.concatenate(per_window, axis=1)
    return u_kn, N_k, per_window


def mbar_leg_dg_kT(u_kn: np.ndarray, N_k: np.ndarray) -> tuple[float, float]:
    """MBAR free energy difference between end states (kT units) and its stat error."""
    mbar = MBAR(u_kn, N_k)
    res = mbar.compute_free_energy_differences()
    return float(res["Delta_f"][0, -1]), float(res["dDelta_f"][0, -1])


def leg_hysteresis_kT(per_window: list[np.ndarray]) -> float:
    """|forward - reverse| Zwanzig estimate across the ladder (kT), a convergence proxy."""
    n_states = len(per_window)
    fwd = rev = 0.0
    for k in range(n_states - 1):
        # forward: samples from state k, d = u_{k+1} - u_k
        uk = per_window[k]
        d_fwd = uk[k + 1] - uk[k]
        fwd += -_logmeanexp(-d_fwd)
        # reverse: samples from state k+1
        ukp1 = per_window[k + 1]
        d_rev = ukp1[k + 1] - ukp1[k]
        rev += _logmeanexp(d_rev)
    return abs(fwd - rev)


def analyze_variant(cfg: dict, variant: str, out_path: str | Path,
                    fep_dir: str | Path | None = None) -> dict:
    """Compute ΔΔG (kcal/mol) + uncertainty + convergence for one variant; write JSON."""
    fep_dir = Path(fep_dir) if fep_dir else ROOT / "results" / "fep"
    fcfg = cfg["fep"]
    n_states = fcfg["lambda_windows"]
    n_reps = fcfg["replicates"]
    legs = fcfg["legs"]
    kT = _KB_KCAL * fcfg["temperature_K"]
    max_closure = fcfg["convergence"]["max_cycle_closure_kcal"]

    per_rep_ddg, per_rep_stat_var, hysteresis_kcal = [], [], []
    for rep in range(n_reps):
        leg_dg, leg_var = {}, {}
        for leg in legs:
            u_kn, N_k, per_window = load_leg_replicate(fep_dir, variant, leg, rep, n_states)
            dg_kT, ddg_kT = mbar_leg_dg_kT(u_kn, N_k)
            leg_dg[leg] = dg_kT * kT
            leg_var[leg] = (ddg_kT * kT) ** 2
            hysteresis_kcal.append(leg_hysteresis_kT(per_window) * kT)
        # ΔΔG = folded - unfolded (positive = destabilizing)
        per_rep_ddg.append(leg_dg["folded"] - leg_dg["unfolded"])
        per_rep_stat_var.append(leg_var["folded"] + leg_var["unfolded"])

    ddg = float(np.mean(per_rep_ddg))
    stat_err = float(np.sqrt(np.sum(per_rep_stat_var)) / n_reps)         # propagated MBAR error
    rep_sem = float(np.std(per_rep_ddg, ddof=1) / np.sqrt(n_reps)) if n_reps > 1 else 0.0
    ddg_err = max(stat_err, rep_sem)                                     # conservative
    cycle_closure = float(np.max(hysteresis_kcal))
    converged = bool(cycle_closure <= max_closure)

    result = {
        "variant": variant,
        "ddg": ddg,
        "ddg_err": ddg_err,
        "cycle_closure_kcal": cycle_closure,
        "converged": converged,
        "n_replicates": n_reps,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="MBAR analysis -> per-variant ddG (Stage 3).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    r = analyze_variant(cfg, args.variant, args.out)
    flag = "" if r["converged"] else "  [NOT CONVERGED]"
    print(f"{args.variant}: ddG = {r['ddg']:.2f} +/- {r['ddg_err']:.2f} kcal/mol "
          f"(cycle closure {r['cycle_closure_kcal']:.2f}){flag} -> {args.out}")


if __name__ == "__main__":
    main()
