"""Stage 3 -- run ONE alchemical FEP window: (variant, leg, window, replicate).

One invocation == one SGE array task == one GPU (README §4 Stage 3). It writes the
reduced potentials for this window to an ``.npz`` that :mod:`src.fep.analyze`
aggregates with MBAR.

The real engine is OpenMM + Perses (GPU) and is isolated behind a backend, so the
pipeline runs locally in **mock mode** (``FEP_MOCK=1``) using synthetic
harmonic-oscillator samples of *known* free energy -- enough to exercise the window
IO, the checkpoint/output schema, and the MBAR analyzer without a GPU. Mock output is
never a result, and that is ENFORCED rather than documented: every window records a
``provenance`` string, :mod:`src.fep.analyze` refuses to mix provenances and stamps the
per-variant ΔΔG with it, and :mod:`src.analysis.validate` refuses to gate on anything
whose provenance is not the configured ``fep.framework``.

NPZ schema (one window):
    lambda_index : int                          state k this window sampled
    n_states     : int                          lambda windows in this leg
    u_kn_window  : (n_states, n_samples) float  reduced potential (U/kT) of this
                                                window's samples at every state
    provenance   : str                          "mock" or the fep.framework that ran it
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

from src.prep.build import load_config
from src.seeds import stable_seed

ROOT = Path(__file__).resolve().parents[2]
_MOCK_SAMPLES = 300
MOCK_PROVENANCE = "mock"


def _stable_target_kT(variant: str, leg: str) -> float:
    """Deterministic per-(variant, leg) free energy (kT) for the mock oscillator ladder.

    Seeded via :func:`src.seeds.stable_seed`, NOT ``hash``: every window of a leg runs in
    its own process, so a per-process seed would give each window a different oscillator
    ladder and MBAR would silently combine incompatible states.
    """
    rng = np.random.default_rng(stable_seed(variant, leg, "dG"))
    return float(rng.uniform(-3.0, 12.0))  # arbitrary but reproducible, kT units


def mock_window(variant: str, leg: str, window: int, rep: int, n_states: int) -> dict:
    """Synthetic 1D harmonic-oscillator samples for one window (analytic free energy).

    State k has U_k(x) = 0.5*K_k*x^2 with K_k = exp(2*dG*k/(n_states-1)); then
    F_k - F_0 = dG*k/(n_states-1), so the ladder's total free energy is exactly ``dG``.
    Samples of state k are drawn as x ~ N(0, 1/sqrt(K_k)); reduced potentials are
    evaluated at every state for MBAR.
    """
    dG = _stable_target_kT(variant, leg)
    k_idx = np.arange(n_states)
    K = np.exp(2.0 * dG * k_idx / (n_states - 1))  # spring constants, K_0 = 1

    rng = np.random.default_rng(stable_seed(variant, leg, window, rep))
    x = rng.normal(0.0, 1.0 / np.sqrt(K[window]), size=_MOCK_SAMPLES)
    # u_kn_window[l, n] = 0.5 * K_l * x_n^2
    u_kn_window = 0.5 * K[:, None] * (x[None, :] ** 2)
    return {"lambda_index": window, "n_states": n_states, "u_kn_window": u_kn_window,
            "provenance": MOCK_PROVENANCE}


def _perses_window(cfg, variant, leg, window, rep, smoke=False):  # pragma: no cover
    """OpenMM+Perses window -- :mod:`src.fep.perses_engine` (SCC-only imports)."""
    from src.fep.perses_engine import run_perses_window
    return run_perses_window(cfg, variant, leg, window, rep, smoke=smoke)


def _pmx_window(cfg, variant, leg, window, rep, smoke=False):  # pragma: no cover
    """GROMACS+pmx window -- :mod:`src.fep.pmx_engine` (SCC-only imports)."""
    from src.fep.pmx_engine import run_pmx_window
    return run_pmx_window(cfg, variant, leg, window, rep, smoke=smoke)


# fep.framework -> engine. Explicit map so an unknown framework fails loudly rather than
# silently running the wrong engine (CLAUDE.md #4).
_ENGINES = {
    "gromacs_pmx": _pmx_window,
    "openmm_perses": _perses_window,
}


def _real_window(cfg, variant, leg, window, rep, smoke=False):  # pragma: no cover
    framework = cfg["fep"]["framework"]
    try:
        engine = _ENGINES[framework]
    except KeyError:
        raise KeyError(
            f"fep.framework={framework!r} has no engine. Known: {sorted(_ENGINES)}."
        ) from None
    return engine(cfg, variant, leg, window, rep, smoke=smoke)


def run_window(cfg: dict, variant: str, leg: str, window: int, rep: int,
               out_path: str | Path, smoke: bool = False) -> Path:
    """Produce one window's reduced potentials and write them to ``out_path`` (.npz)."""
    n_states = cfg["fep"]["lambda_windows"]
    if os.environ.get("FEP_MOCK"):
        data = mock_window(variant, leg, window, rep, n_states)
    else:
        data = _real_window(cfg, variant, leg, window, rep, smoke=smoke)
    # A window without provenance cannot be told apart from a hand-written file.
    data.setdefault("provenance", MOCK_PROVENANCE if os.environ.get("FEP_MOCK")
                    else cfg["fep"]["framework"])

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, **data)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one alchemical FEP window (Stage 3).")
    parser.add_argument("--variant", required=True)
    parser.add_argument("--leg", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--smoke", action="store_true",
                        help="Tiny step/frame counts to shake out the GPU/Perses path fast "
                             "(real window; not a result). Ignored under FEP_MOCK.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out = run_window(cfg, args.variant, args.leg, args.window, args.rep, args.out,
                     smoke=args.smoke)
    print(f"Wrote window {args.variant}/{args.leg}/w{args.window}_r{args.rep} -> {out}")


if __name__ == "__main__":
    main()
