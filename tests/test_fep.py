"""Stage 3 FEP: window IO (mock) and the MBAR analyzer (real pymbar).

The analyzer is validated against synthetic harmonic-oscillator ladders of KNOWN
free energy, so MBAR must recover the injected ΔΔG. No GPU/Perses required.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from src.fep.analyze import _KB_KCAL, analyze_variant, leg_hysteresis_kT
from src.fep.window import run_window

FEP_CFG = {
    "fep": {
        "lambda_windows": 12,
        "replicates": 3,
        "legs": ["folded", "unfolded"],
        "temperature_K": 298.15,
        "convergence": {"max_cycle_closure_kcal": 1.0},
    }
}
_KT = _KB_KCAL * 298.15


def _write_oscillator_leg(fep_dir: Path, variant: str, leg: str, dG_kT: float,
                          n_states: int, n_reps: int, n_samples: int = 400) -> None:
    """Write a leg's windows as harmonic oscillators with total free energy dG_kT."""
    K = np.exp(2.0 * dG_kT * np.arange(n_states) / (n_states - 1))  # K_0 = 1
    for rep in range(n_reps):
        rng = np.random.default_rng(abs(hash((variant, leg, rep))) % (2**32))
        for w in range(n_states):
            x = rng.normal(0.0, 1.0 / np.sqrt(K[w]), size=n_samples)
            u_kn_window = 0.5 * K[:, None] * (x[None, :] ** 2)
            d = fep_dir / variant / leg
            d.mkdir(parents=True, exist_ok=True)
            np.savez(d / f"w{w}_r{rep}.npz",
                     lambda_index=w, n_states=n_states, u_kn_window=u_kn_window)


def test_run_window_mock_schema(tmp_path, monkeypatch):
    monkeypatch.setenv("FEP_MOCK", "1")
    out = tmp_path / "w0_r0.npz"
    run_window(FEP_CFG, "A4V", "folded", 3, 0, out)
    npz = np.load(out)
    assert int(npz["n_states"]) == 12
    assert int(npz["lambda_index"]) == 3
    assert npz["u_kn_window"].shape[0] == 12
    assert np.all(np.isfinite(npz["u_kn_window"]))


def test_mbar_recovers_known_ddg(tmp_path):
    # ΔΔG = (folded - unfolded) = (3.0 - 1.0) kT
    fep_dir = tmp_path / "fep"
    _write_oscillator_leg(fep_dir, "TEST", "folded", 3.0, 12, 3)
    _write_oscillator_leg(fep_dir, "TEST", "unfolded", 1.0, 12, 3)

    result = analyze_variant(FEP_CFG, "TEST", tmp_path / "ddg.json", fep_dir=fep_dir)
    expected = (3.0 - 1.0) * _KT
    assert result["ddg"] == pytest.approx(expected, abs=0.2)  # MBAR recovers injected ddG
    assert result["ddg_err"] > 0.0
    assert result["converged"] is True
    assert result["cycle_closure_kcal"] < 1.0


def test_analyze_writes_gate_schema(tmp_path):
    fep_dir = tmp_path / "fep"
    _write_oscillator_leg(fep_dir, "TEST", "folded", 2.0, 12, 3)
    _write_oscillator_leg(fep_dir, "TEST", "unfolded", 2.0, 12, 3)
    out = tmp_path / "ddg.json"
    analyze_variant(FEP_CFG, "TEST", out, fep_dir=fep_dir)
    saved = json.loads(out.read_text())
    assert set(saved) >= {"variant", "ddg", "ddg_err", "cycle_closure_kcal", "converged"}
    # equal legs -> ΔΔG ~ 0
    assert saved["ddg"] == pytest.approx(0.0, abs=0.2)


def test_hysteresis_small_for_good_overlap():
    # a single well-sampled ladder should have near-zero forward/reverse hysteresis
    n_states = 12
    K = np.exp(2.0 * 2.0 * np.arange(n_states) / (n_states - 1))
    rng = np.random.default_rng(0)
    per_window = []
    for w in range(n_states):
        x = rng.normal(0.0, 1.0 / np.sqrt(K[w]), size=2000)
        per_window.append(0.5 * K[:, None] * (x[None, :] ** 2))
    assert leg_hysteresis_kT(per_window) < 0.5  # kT
