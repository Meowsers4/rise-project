"""Stage 3 (analysis) -- windows -> per-variant ΔΔG via MBAR, with uncertainty.

Aggregates the per-window reduced potentials (from :mod:`src.fep.window`) over both
legs and all replicates, runs MBAR per (leg, replicate), and forms

    ΔΔG_folding = ΔG(folded leg) - ΔG(unfolded leg)      (positive = destabilizing)

averaged over replicates. A ΔΔG is reported only WITH an uncertainty and a
convergence check (CLAUDE.md #5): the forward/reverse (Zwanzig) hysteresis of each
leg is the available per-variant convergence proxy, and ``converged`` is False if it
exceeds ``fep.convergence.max_cycle_closure_kcal``.

Provenance is carried, not assumed: every window states which engine produced it, all
windows of a variant must agree, and the value is stamped onto the ΔΔG so the gate can
refuse anything that is not a production run. A window with no provenance is rejected --
that is exactly what a hand-written or hand-edited file looks like.

Emits ``results/fep/<variant>/ddg.json`` in the schema the validation gate consumes
(see :data:`src.analysis.validate.FEP_DDG_SCHEMA`).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from src.prep.build import load_config

ROOT = Path(__file__).resolve().parents[2]
_KB_KCAL = 0.0019872041  # Boltzmann constant, kcal/mol/K


def _logmeanexp(a: np.ndarray) -> float:
    """log(mean(exp(a))) with the max-shift trick (avoids overflow)."""
    m = np.max(a)
    return float(m + np.log(np.mean(np.exp(a - m))))


_MIN_SAMPLES_TO_DECORRELATE = 10      # below this, timeseries analysis is noise; keep as-is

# Windows written before the .mdp fingerprint existed. Treated as its own protocol so a
# pre-fingerprint set still analyses alone but can never be silently mixed with a later one.
_UNLABELLED_PROTOCOL = "unlabelled"


def decorrelate_window(u: np.ndarray, lambda_index: int) -> np.ndarray:
    """Trim equilibration and subsample ``u`` to statistically independent frames.

    MBAR's uncertainty assumes independent samples. Frames saved every 20 ps are not:
    the system has not had time to forget where it was, so N frames carry the
    information of roughly N/g independent ones and MBAR's error comes out about
    ``sqrt(g)`` times too small. ``pymbar.timeseries`` estimates the equilibration point
    and the correlation time g from the window's OWN reduced potential (the observable
    it actually sampled), and we keep only the independent frames.

    Args:
        u: ``(n_states, n_samples)`` reduced potentials for one window.
        lambda_index: Which row of ``u`` is the state this window sampled.

    Returns:
        ``(n_states, n_independent)`` -- the same rows, decorrelated columns.
    """
    if u.shape[1] < _MIN_SAMPLES_TO_DECORRELATE:
        return u
    from pymbar import timeseries

    u_self = u[lambda_index]                       # potential of the state actually sampled
    # u_self is the TOTAL reduced potential, dominated by ~30k solvent atoms whose
    # autocorrelation is sub-ps. On its own it reports g ~ 1 at any sane frame spacing --
    # which looks like "our samples are independent" while saying nothing about the
    # alchemical coordinate, the slow mode that actually limits this calculation. Take the
    # neighbour energy difference too (what MBAR consumes) and keep the LARGER g, so the
    # thinning can only ever be more conservative, never less.
    nbr = lambda_index + 1 if lambda_index + 1 < u.shape[0] else lambda_index - 1
    du = u[nbr] - u_self if nbr >= 0 else u_self
    t0, g_self, _ = timeseries.detect_equilibration(u_self)
    t0_du, g_du, _ = timeseries.detect_equilibration(du)
    t0, g = max(t0, t0_du), max(g_self, g_du)
    keep = timeseries.subsample_correlated_data(u_self[t0:], g=g)
    return u[:, t0:][:, keep]


def load_leg_replicate(fep_dir: Path, variant: str, leg: str, rep: int, n_states: int,
                       decorrelate: bool = True
                       ) -> tuple[np.ndarray, np.ndarray, list[np.ndarray], int]:
    """Aggregate a leg's windows for one replicate into MBAR inputs.

    Returns (u_kn, N_k, per_window_u, n_raw) where u_kn is (n_states, total_samples),
    per_window_u[k] is window k's (n_states, n_k) block (used for hysteresis), and n_raw
    is the sample count before decorrelation (so the thinning is visible in the output).

    Raises:
        ValueError: If a window's ladder does not match the configured state count.
    """
    per_window = []
    n_raw = 0
    for w in range(n_states):
        npz = np.load(fep_dir / variant / leg / f"w{w}_r{rep}.npz")
        u = npz["u_kn_window"]
        if u.shape[0] != n_states:
            raise ValueError(
                f"{variant}/{leg}/w{w}_r{rep}: u_kn_window has {u.shape[0]} states but "
                f"config says {n_states} -- these windows are not the same ladder."
            )
        n_raw += u.shape[1]
        if decorrelate:
            u = decorrelate_window(u, int(npz["lambda_index"]))
        per_window.append(u)
    N_k = np.array([u.shape[1] for u in per_window])
    u_kn = np.concatenate(per_window, axis=1)
    return u_kn, N_k, per_window, n_raw


def collect_provenance(fep_dir: Path, variant: str, legs: list[str], n_reps: int,
                       n_states: int) -> tuple[set[str], set[str]]:
    """Engine tags across every window of ``variant``, read without loading the samples.

    Runs before any MBAR work so a provenance problem is reported immediately rather than
    after minutes of computation on numbers that were never going to be usable.

    Returns:
        ``(provenances, protocols)``. Provenance says *which engine*; the protocol
        fingerprint says *which settings*. Windows can share an engine and still be
        unmixable -- e.g. half a variant run before ``sc-coul`` was enabled and half
        after. Returned rather than collected through an out-parameter so a caller
        cannot drop the protocol set and silently disable the mixing guard.

    Raises:
        ValueError: If any window carries no provenance -- which is what a hand-written
            or pre-provenance file looks like, and is never trustworthy.
    """
    provenances: set[str] = set()
    protocols: set[str] = set()
    for leg in legs:
        for rep in range(n_reps):
            for w in range(n_states):
                path = fep_dir / variant / leg / f"w{w}_r{rep}.npz"
                with np.load(path) as npz:      # lazy: only the provenance array is read
                    if "provenance" not in npz:
                        raise ValueError(
                            f"{path}: no provenance recorded. A window written before "
                            "provenance tracking (or written by hand) cannot be "
                            "distinguished from a real run; re-run it with src.fep.window."
                        )
                    provenances.add(str(npz["provenance"]))
                    protocols.add(str(npz["protocol"]) if "protocol" in npz
                                  else _UNLABELLED_PROTOCOL)
    return provenances, protocols


def _solver_notes(chatter: str) -> set[str]:
    """Solver chatter worth recording, with import-time banners stripped.

    pymbar's JAX backend prints a six-line boxed "JAX 64-bit mode is now on!" banner to
    stderr the first time it is imported. It carries no information about the solve, and
    it dominated ``solver_notes`` on the F64A run (five unique lines, the whole field).
    Nothing was lost -- the notes are an unbounded set, so a real fallback warning would
    still have appeared alongside it -- but a diagnostics field that is pure boilerplate
    does not get read, which defeats the reason 8800b73 started capturing stderr at all.

    Drops banner box-drawing (lines that are only ``*`` and spaces) and any line inside
    such a box, keeping anything pymbar says about the solve itself.
    """
    notes: set[str] = set()
    for raw in chatter.splitlines():
        line = raw.strip()
        if not line:
            continue
        # Banner lines are wrapped in '*' or consist only of '*' padding.
        if line.startswith("*") and line.endswith("*"):
            continue
        if set(line) <= {"*", " "}:
            continue
        notes.add(line)
    return notes


def solve_leg_mbar(u_kn: np.ndarray, N_k: np.ndarray) -> dict:
    """Solve MBAR ONCE for a leg and return everything downstream needs.

    Previously ``mbar_leg_dg_kT`` and ``mbar_leg_overlap`` each constructed their own
    ``MBAR``, so every (leg, replicate) was solved twice. Harmless at 151 samples/window;
    at 3001 it doubles an already expensive step, and it also allowed the free energy and
    the overlap matrix to come from two independently-converged solves.

    Captures pymbar's solver chatter. "Failed to reach a solution to within tolerance with
    hybr: trying next method" means the first self-consistent solver gave up and a
    fallback took over -- pymbar prints this to stdout rather than raising, so without
    capturing it a marginal solve is invisible. Recorded, not fatal: the fallbacks are
    legitimate, but a run that needed them should say so in its diagnostics.

    Returns:
        dict with ``dg_kT``, ``ddg_kT``, ``overlap`` (nested lists), ``min_adjacent``
        and ``solver_notes``.

    Raises:
        FloatingPointError: If MBAR returns a non-finite free energy or uncertainty --
            an unconverged solve that would otherwise become a confident-looking number.
    """
    import contextlib
    import io

    from pymbar import MBAR  # lazy: keeps the pure helpers importable without pymbar

    # BOTH streams: pymbar prints "Failed to reach a solution to within tolerance with
    # hybr: trying next method" to STDERR, so capturing stdout alone reported an empty
    # solver_notes for a G93A run that had five fallbacks -- instrumentation that says
    # "clean" when it is not is worse than none.
    chatter = io.StringIO()
    with contextlib.redirect_stdout(chatter), contextlib.redirect_stderr(chatter):
        mbar = MBAR(u_kn, N_k)
        res = mbar.compute_free_energy_differences()
        overlap = np.asarray(mbar.compute_overlap()["matrix"], dtype=float)

    dg = float(res["Delta_f"][0, -1])
    ddg = float(res["dDelta_f"][0, -1])
    if not (np.isfinite(dg) and np.isfinite(ddg)):
        raise FloatingPointError(
            f"MBAR returned a non-finite free energy (dg={dg}, ddg={ddg}). The solver did "
            "not converge; refusing to report a number."
        )

    adjacent = [float(overlap[k, k + 1]) for k in range(overlap.shape[0] - 1)]
    notes = sorted(_solver_notes(chatter.getvalue()))
    return {
        "dg_kT": dg,
        "ddg_kT": ddg,
        "overlap": overlap.tolist(),
        "min_adjacent": min(adjacent) if adjacent else float("nan"),
        "solver_notes": notes,
    }


def mbar_leg_dg_kT(u_kn: np.ndarray, N_k: np.ndarray) -> tuple[float, float]:
    """Free energy difference between end states (kT) and its stat error.

    Thin wrapper on :func:`solve_leg_mbar`, kept because it is the documented unit of the
    analyzer and is what the oscillator tests exercise.
    """
    r = solve_leg_mbar(u_kn, N_k)
    return r["dg_kT"], r["ddg_kT"]


def mbar_leg_overlap(u_kn: np.ndarray, N_k: np.ndarray) -> tuple[list[list[float]], float]:
    """MBAR phase-space overlap matrix and the worst adjacent-window overlap.

    The nearest-neighbour diagonal decides whether a lambda ladder is dense enough: if
    adjacent windows barely overlap, MBAR is extrapolating and its error bar understates
    the truth. Carries claim C2, which cannot be argued retroactively.
    """
    r = solve_leg_mbar(u_kn, N_k)
    return r["overlap"], r["min_adjacent"]


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


def _check_single_protocol(variant: str, protocols: set[str]) -> None:
    """Raise unless every window was produced by the same .mdp settings.

    The engine tag cannot see this: `sc-coul` off vs on, or a changed `nstdhdl`, are the
    same engine producing incompatible samples. Combining them in one MBAR estimate is
    silent and wrong. Windows written before fingerprinting count as their own protocol,
    so a pre-change set still analyses on its own but can never be mixed with a later one.

    Raises:
        ValueError: If the windows disagree on the mdp fingerprint.
    """
    if len(protocols) > 1:
        raise ValueError(
            f"{variant}: windows were run under {len(protocols)} different protocols "
            f"{sorted(protocols)}. A ΔΔG may not mix .mdp settings (sc-coul, nstdhdl, "
            "sampling length ...) any more than it may mix engines. Re-run the variant "
            "so every window shares one protocol."
        )


def _check_single_provenance(variant: str, provenances: set[str]) -> None:
    """Raise unless every window seen so far came from the same engine.

    Raises:
        ValueError: If the windows mix engines (e.g. a mock window quietly standing in
            for a GPU window that failed).
    """
    if len(provenances) > 1:
        raise ValueError(
            f"{variant}: windows disagree on provenance {sorted(provenances)}. A ΔΔG may "
            "not mix engines -- re-run the odd windows so the variant comes from one."
        )


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

    # Establish provenance for the whole variant BEFORE computing anything with it.
    provenances, protocols = collect_provenance(fep_dir, variant, legs, n_reps, n_states)
    _check_single_provenance(variant, provenances)
    _check_single_protocol(variant, protocols)
    provenance = provenances.pop()
    protocol = protocols.pop() if protocols else _UNLABELLED_PROTOCOL

    decorrelate = bool(fcfg.get("decorrelate", True))
    per_rep_ddg, per_rep_stat_var, hysteresis_kcal = [], [], []
    diagnostics: list[dict] = []
    solver_notes: set[str] = set()
    n_raw_total = n_used_total = 0
    for rep in range(n_reps):
        leg_dg, leg_var = {}, {}
        for leg in legs:
            u_kn, N_k, per_window, n_raw = load_leg_replicate(
                fep_dir, variant, leg, rep, n_states, decorrelate=decorrelate)
            n_raw_total += n_raw
            n_used_total += int(N_k.sum())
            solved = solve_leg_mbar(u_kn, N_k)          # one MBAR per (leg, replicate)
            leg_dg[leg] = solved["dg_kT"] * kT
            leg_var[leg] = (solved["ddg_kT"] * kT) ** 2
            leg_hyst = leg_hysteresis_kT(per_window) * kT
            hysteresis_kcal.append(leg_hyst)
            if solved["solver_notes"]:
                solver_notes.update(solved["solver_notes"])
            diagnostics.append({
                "leg": leg, "replicate": rep,
                "dg_kcal": leg_dg[leg], "dg_err_kcal": float(solved["ddg_kT"] * kT),
                "hysteresis_kcal": leg_hyst,
                "min_adjacent_overlap": solved["min_adjacent"],
                "samples_per_window": [int(n) for n in N_k],
                "solver_notes": solved["solver_notes"],
                "overlap_matrix": solved["overlap"],
            })
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
        "provenance": provenance,
        "protocol": protocol,
        "decorrelated": decorrelate,
        "n_samples_raw": n_raw_total,
        "n_samples_independent": n_used_total,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2))

    # Claim C2 depends on these being persisted at run time, not reconstructed later
    # (README §2.3, §5). One file per variant, alongside the ddG. Sited relative to
    # fep_dir -- NOT to ROOT -- so a test pointing fep_dir at a tmp path does not write
    # into the real results tree.
    conv_dir = fep_dir.parent / "convergence"
    conv_dir.mkdir(parents=True, exist_ok=True)
    (conv_dir / f"{variant}.json").write_text(json.dumps({
        "variant": variant,
        "provenance": provenance,
        "protocol": protocol,
        "ddg": ddg,
        "ddg_err": ddg_err,
        "per_replicate_ddg": [float(x) for x in per_rep_ddg],
        "replicate_spread_kcal": float(np.ptp(per_rep_ddg)) if n_reps > 1 else 0.0,
        "replicate_sem_kcal": rep_sem,
        "cycle_closure_kcal": cycle_closure,
        "max_cycle_closure_kcal": max_closure,
        "converged": converged,
        "min_adjacent_overlap": min(d["min_adjacent_overlap"] for d in diagnostics),
        "solver_notes": sorted(solver_notes),
        "legs": diagnostics,
    }, indent=2))
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
    if r["provenance"] != cfg["fep"]["framework"]:
        flag += f"  [{r['provenance'].upper()} -- NOT A RESULT, the gate will reject this]"
    print(f"{args.variant}: ddG = {r['ddg']:.2f} +/- {r['ddg_err']:.2f} kcal/mol "
          f"(cycle closure {r['cycle_closure_kcal']:.2f}){flag} -> {args.out}")
    if r["decorrelated"]:
        print(f"  independent samples: {r['n_samples_independent']} of "
              f"{r['n_samples_raw']} raw frames")


if __name__ == "__main__":
    main()
