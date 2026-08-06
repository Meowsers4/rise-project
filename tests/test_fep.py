"""Stage 3 FEP: window IO (mock) and the MBAR analyzer (real pymbar).

The analyzer is validated against synthetic harmonic-oscillator ladders of KNOWN
free energy, so MBAR must recover the injected ΔΔG. No GPU/Perses required.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from src.fep.analyze import _KB_KCAL, analyze_variant, leg_hysteresis_kT
from src.fep.window import run_window
from src.seeds import stable_seed

FEP_CFG = {
    "fep": {
        "framework": "openmm_perses",
        "lambda_windows": 12,
        "replicates": 3,
        "legs": ["folded", "unfolded"],
        "temperature_K": 298.15,
        "convergence": {"max_cycle_closure_kcal": 1.0},
    }
}
_KT = _KB_KCAL * 298.15


def _write_oscillator_leg(fep_dir: Path, variant: str, leg: str, dG_kT: float,
                          n_states: int, n_reps: int, n_samples: int = 400,
                          provenance: str = "openmm_perses") -> None:
    """Write a leg's windows as harmonic oscillators with total free energy dG_kT."""
    K = np.exp(2.0 * dG_kT * np.arange(n_states) / (n_states - 1))  # K_0 = 1
    for rep in range(n_reps):
        rng = np.random.default_rng(stable_seed(variant, leg, rep))
        for w in range(n_states):
            x = rng.normal(0.0, 1.0 / np.sqrt(K[w]), size=n_samples)
            u_kn_window = 0.5 * K[:, None] * (x[None, :] ** 2)
            d = fep_dir / variant / leg
            d.mkdir(parents=True, exist_ok=True)
            np.savez(d / f"w{w}_r{rep}.npz",
                     lambda_index=w, n_states=n_states, u_kn_window=u_kn_window,
                     provenance=provenance)


def test_run_window_mock_schema(tmp_path, monkeypatch):
    monkeypatch.setenv("FEP_MOCK", "1")
    out = tmp_path / "w0_r0.npz"
    run_window(FEP_CFG, "A4V", "folded", 3, 0, out)
    npz = np.load(out)
    assert int(npz["n_states"]) == 12
    assert int(npz["lambda_index"]) == 3
    assert npz["u_kn_window"].shape[0] == 12
    assert np.all(np.isfinite(npz["u_kn_window"]))
    assert str(npz["provenance"]) == "mock"   # mock must be self-identifying


def test_real_window_requires_validation_gate(tmp_path, monkeypatch):
    monkeypatch.delenv("FEP_MOCK", raising=False)
    cfg = {
        **FEP_CFG,
        "validation": {
            "gate_subset": ["A4V"],
            "outputs": {"gate_report": str(tmp_path / "missing-gate.json")},
        },
    }
    with pytest.raises(RuntimeError, match="validation gate"):
        run_window(cfg, "G93A", "folded", 0, 0, tmp_path / "window.npz")


def test_mock_window_is_identical_across_processes():
    """Regression: hash() is per-process randomized, so seeding from it gave every SGE
    array task a different oscillator ladder and MBAR silently blended them."""
    root = Path(__file__).resolve().parents[1]
    code = ("from src.fep.window import _stable_target_kT, mock_window;"
            "d=mock_window('A4V','folded',3,0,12);"
            "print(_stable_target_kT('A4V','folded'), d['u_kn_window'].sum())")
    outs = set()
    for hashseed in ("0", "1", "random"):
        env = {**os.environ, "PYTHONHASHSEED": hashseed}
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                           cwd=root, env=env, check=True)
        outs.add(r.stdout.strip())
    assert len(outs) == 1, f"mock window is not process-stable: {outs}"


def test_analyze_rejects_mixed_provenance(tmp_path):
    fep_dir = tmp_path / "fep"
    _write_oscillator_leg(fep_dir, "TEST", "folded", 2.0, 12, 1)
    _write_oscillator_leg(fep_dir, "TEST", "unfolded", 1.0, 12, 1, provenance="mock")
    cfg = {"fep": {**FEP_CFG["fep"], "replicates": 1}}
    with pytest.raises(ValueError, match="disagree on provenance"):
        analyze_variant(cfg, "TEST", tmp_path / "ddg.json", fep_dir=fep_dir)


def test_analyze_rejects_window_without_provenance(tmp_path):
    """An unlabelled window is indistinguishable from a hand-written one."""
    fep_dir = tmp_path / "fep"
    _write_oscillator_leg(fep_dir, "TEST", "folded", 2.0, 12, 1)
    _write_oscillator_leg(fep_dir, "TEST", "unfolded", 1.0, 12, 1)
    stripped = fep_dir / "TEST" / "folded" / "w0_r0.npz"
    d = dict(np.load(stripped))
    d.pop("provenance")
    np.savez(stripped, **d)
    cfg = {"fep": {**FEP_CFG["fep"], "replicates": 1}}
    with pytest.raises(ValueError, match="no provenance recorded"):
        analyze_variant(cfg, "TEST", tmp_path / "ddg.json", fep_dir=fep_dir)


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
    assert set(saved) >= {"variant", "ddg", "ddg_err", "cycle_closure_kcal", "converged",
                          "provenance"}
    assert saved["provenance"] == "openmm_perses"   # stamped through from the windows
    # equal legs -> ΔΔG ~ 0
    assert saved["ddg"] == pytest.approx(0.0, abs=0.2)


def test_decorrelation_thins_correlated_samples(tmp_path):
    """Correlated frames must be thinned; i.i.d. frames should survive nearly intact."""
    from src.fep.analyze import decorrelate_window

    n_states, n = 6, 600
    rng = np.random.default_rng(0)
    # strongly autocorrelated series (AR(1), phi=0.95) -> should thin hard
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = 0.95 * x[i - 1] + rng.normal(0, 1)
    u_corr = np.tile(x, (n_states, 1))
    kept_corr = decorrelate_window(u_corr, 0).shape[1]

    u_iid = np.tile(rng.normal(0, 1, n), (n_states, 1))
    kept_iid = decorrelate_window(u_iid, 0).shape[1]

    assert kept_corr < n / 4, f"correlated series barely thinned ({kept_corr}/{n})"
    assert kept_iid > n / 2, f"i.i.d. series over-thinned ({kept_iid}/{n})"


def test_decorrelation_recorded_in_output(tmp_path):
    fep_dir = tmp_path / "fep"
    _write_oscillator_leg(fep_dir, "TEST", "folded", 2.0, 12, 1)
    _write_oscillator_leg(fep_dir, "TEST", "unfolded", 2.0, 12, 1)
    cfg = {"fep": {**FEP_CFG["fep"], "replicates": 1, "decorrelate": True}}
    r = analyze_variant(cfg, "TEST", tmp_path / "ddg.json", fep_dir=fep_dir)
    assert r["decorrelated"] is True
    assert 0 < r["n_samples_independent"] <= r["n_samples_raw"]


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


def _fake_pmx(monkeypatch, tmp_path, layout: dict[str, list[str]]):
    """Install a stub ``pmx`` package whose data dir has the given <dir>/<ff>.ff layout."""
    import types

    pkg = tmp_path / "pmx"
    for dirname, ffs in layout.items():
        for ff in ffs:
            (pkg / "data" / dirname / f"{ff}.ff").mkdir(parents=True, exist_ok=True)
    module = types.ModuleType("pmx")
    module.__file__ = str(pkg / "__init__.py")
    monkeypatch.setitem(sys.modules, "pmx", module)


def _ff_cfg(name="amber99sb-star-ildn-mut"):
    return {"fep": {"pmx_ff_dir": "auto", "pmx_forcefield": name}}


def test_pmx_ff_dir_prefers_mutff45_on_old_layout(monkeypatch, tmp_path):
    """pmx 2.0 ships both; mutff45 holds the real set and must win over legacy mutff."""
    from src.fep.pmx_engine import pmx_ff_dir

    _fake_pmx(monkeypatch, tmp_path, {
        "mutff45": ["amber99sb-star-ildn-mut"],
        "mutff": ["amber99sb-star-ildn-mut"],
    })
    assert pmx_ff_dir(_ff_cfg()).name == "mutff45"


def test_pmx_ff_dir_finds_mutff_on_current_develop(monkeypatch, tmp_path):
    """Current develop dropped mutff45; the force field lives in mutff."""
    from src.fep.pmx_engine import pmx_ff_dir

    _fake_pmx(monkeypatch, tmp_path, {"mutff": ["amber99sb-star-ildn-mut"]})
    assert pmx_ff_dir(_ff_cfg()).name == "mutff"


def test_pmx_ff_dir_skips_dir_lacking_the_forcefield(monkeypatch, tmp_path):
    """A mutff45 that does not carry the configured ff must not shadow the one that does."""
    from src.fep.pmx_engine import pmx_ff_dir

    _fake_pmx(monkeypatch, tmp_path, {
        "mutff45": ["charmm36mut"],
        "mutff": ["amber99sb-star-ildn-mut"],
    })
    assert pmx_ff_dir(_ff_cfg()).name == "mutff"


def test_pmx_ff_dir_explicit_config_wins(monkeypatch, tmp_path):
    from src.fep.pmx_engine import pmx_ff_dir

    _fake_pmx(monkeypatch, tmp_path, {"mutff": ["amber99sb-star-ildn-mut"]})
    cfg = {"fep": {"pmx_ff_dir": "/opt/custom/ff", "pmx_forcefield": "x"}}
    assert pmx_ff_dir(cfg) == Path("/opt/custom/ff")
