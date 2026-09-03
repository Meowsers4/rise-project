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
ROOT = Path(__file__).resolve().parents[1]


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


def test_analyze_refuses_to_mix_protocols_end_to_end(tmp_path):
    """The guard must be WIRED, not merely present.

    An earlier version collected protocols through a mutable out-parameter; deleting the
    argument at the call site disabled the whole check and every test still passed. This
    goes through analyze_variant with real files on disk, so the collection path itself
    is exercised.
    """
    from src.fep.analyze import analyze_variant

    cfg = {**FEP_CFG, "fep": {**FEP_CFG["fep"], "lambda_windows": 4, "replicates": 1,
                              "decorrelate": False}}
    fep_dir = tmp_path / "fep"
    for leg in ("folded", "unfolded"):
        d = fep_dir / "A4V" / leg
        d.mkdir(parents=True)
        for w in range(4):
            # one window disagrees about the .mdp it was produced with
            proto = "PULLED_MIDWAY" if (leg, w) == ("folded", 2) else "abc123"
            np.savez(d / f"w{w}_r0.npz", lambda_index=w, n_states=4,
                     u_kn_window=np.zeros((4, 50)), provenance="openmm_perses",
                     protocol=proto)
    with pytest.raises(ValueError, match="different protocols"):
        analyze_variant(cfg, "A4V", tmp_path / "ddg.json", fep_dir=fep_dir)
    assert not (tmp_path / "ddg.json").exists()   # and no result is written


def test_analyze_refuses_to_mix_protocols():
    """Same engine, different .mdp settings, is still unmixable -- WITHIN a leg.

    Provenance catches a mock or another engine. It cannot see sc-coul flipped or nstdhdl
    changed halfway through a leg -- which is what `git pull`ing the cluster during a
    running array produces.

    Across legs is different: fep.ns_per_window and fep.equilibration_ns may be per-leg,
    so the folded and unfolded legs hash differently BY DESIGN. Each leg is its own MBAR
    estimate and they never share samples, so that must be allowed.
    """
    from src.fep.analyze import _check_single_protocol

    ok = {"folded": {"abc123"}, "unfolded": {"abc123"}}
    _check_single_protocol("A4V", ok)

    # legs may legitimately differ -- per-leg sampling
    per_leg = {"folded": {"folded9ns"}, "unfolded": {"unfolded3ns"}}
    _check_single_protocol("A4V", per_leg)

    # but within one leg they may not
    with pytest.raises(ValueError, match="A4V/folded"):
        _check_single_protocol("A4V", {"folded": {"abc123", "def456"},
                                       "unfolded": {"abc123"}})
    with pytest.raises(ValueError, match="A4V/unfolded"):
        _check_single_protocol("A4V", {"folded": {"abc123"},
                                       "unfolded": {"abc123", "unlabelled"}})


def test_forcefield_protonation_names_are_not_wrong_residues():
    """pdb2gmx writes HID/HIE/HIP, never HIS; the panel speaks in PDB names.

    Without aliasing, every histidine variant in the panel (H43R, H46R, H71Q, H110Y,
    H120L) dies with a spurious 'refusing to mutate the wrong residue'.
    """
    from src.fep.pmx_engine import _canonical_resname

    for ff_name in ("HID", "HIE", "HIP", "HSD", "HSE", "HSP"):
        assert _canonical_resname(ff_name) == "HIS"
    assert _canonical_resname("CYX") == "CYS"
    assert _canonical_resname("ASH") == "ASP"
    assert _canonical_resname("VAL") == "VAL"      # ordinary names pass through
    assert _canonical_resname(None) is None


def test_window_minimisation_is_lambda_aware():
    """Regression: EM with free-energy OFF minimises at the A state, where the appearing
    atoms are dummies. Real atoms relax on top of them and the window detonates on step 1
    at high lambda (A4V folded w17, all three replicates, 2026-08-06)."""
    from src.fep.pmx_engine import _minim_mdp

    cfg = {"fep": {"lambda_windows": 18, "nonbonded_method": "PME",
                   "nonbonded_cutoff_nm": 1.0, "sc_function": "beutler",
                   "sc_alpha": 0.3, "sc_power": 1, "sc_sigma": 0.25}}

    per_window = _minim_mdp(cfg, 17)
    assert "free-energy              = yes" in per_window
    assert "init-lambda-state        = 17" in per_window
    assert "sc-function              = beutler" in per_window   # soft-core, or lambda=1 blows up

    # The shared system build minimises the physical topology only -- no lambda there.
    assert "free-energy" not in _minim_mdp(cfg)


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


def test_protocol_fingerprint_is_seed_and_window_invariant_but_settings_sensitive():
    """Guards the exclusion list in run_pmx_window's fingerprint.

    If the excluded lines ever stop matching, every replicate hashes differently and
    _check_single_protocol rejects EVERY variant -- after the GPU time is already spent.
    If they exclude too much, a real settings change stops being detected.
    """
    import copy
    import yaml
    from src.fep.pmx_engine import production_mdp, protocol_fingerprint

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))

    def fp(c, window, seed):
        # call the PRODUCTION helper -- an earlier version of this test recomputed the
        # hash inline, so breaking the exclusion list in the engine changed nothing here.
        return protocol_fingerprint(production_mdp(c, window, seed, leg="folded"))

    assert fp(cfg, 0, 1) == fp(cfg, 17, 999)        # replicates/windows must agree

    off = copy.deepcopy(cfg)
    off["fep"]["sc_alpha"] = 0.5                    # a real protocol change must not
    assert fp(off, 0, 1) != fp(cfg, 0, 1)           # hash the same

    coarse = copy.deepcopy(cfg)
    coarse["fep"]["frames_per_window"] = 150
    assert fp(coarse, 0, 1) != fp(cfg, 0, 1)


def test_pre_registered_thresholds_are_not_quietly_lowered():
    """CLAUDE.md rule 2 / HANDOFF standing rule, enforced mechanically.

    These were pre-registered 2026-08-07 before any gate evaluation. Until this test
    existed the suite passed with min_pearson set to 0.10.
    """
    import yaml

    v = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))["validation"]
    f = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))["fep"]
    assert v["min_pearson"] == 0.70
    assert v["max_rmse_kcal"] == 1.5
    assert v["max_median_cycle_closure_kcal"] == 0.75
    assert v["pivot_pearson"] == 0.60
    # internal coherence: the pivot line sits below the pass mark, and the set-level
    # hysteresis bound is tighter than the per-variant cap it complements.
    assert v["pivot_pearson"] < v["min_pearson"]
    assert v["max_median_cycle_closure_kcal"] <= f["convergence"]["max_cycle_closure_kcal"]


def test_disulfide_guard_catches_a_rebuilt_bond(tmp_path):
    """CLAUDE.md: the GROMACS path needs its OWN disulfide guard.

    Until now it piped "n" into pdb2gmx -ss open-loop and never checked. A re-formed
    Cys57-Cys146 bond still yields a finite, plausible ddG, so nothing downstream notices.
    """
    from src.fep.pmx_engine import assert_topology_disulfide_free

    def gro(lines):
        p = tmp_path / f"c{len(list(tmp_path.iterdir()))}.gro"
        p.write_text("t\n%d\n" % len(lines) + "\n".join(lines) + "\n 3 3 3\n")
        return p

    reduced = gro(["   57CYS     SG   1   1.000   1.000   1.000",
                   "   57CYS     HG   2   1.100   1.000   1.000",
                   "  146CYS     SG   3   2.000   2.000   2.000",
                   "  146CYS     HG   4   2.100   2.000   2.000"])
    clean_top = tmp_path / "ok.top"
    clean_top.write_text("[ moleculetype ]\nProtein 3\n; CYS 57\n")
    assert_topology_disulfide_free(clean_top, reduced)      # must not raise

    # 1. bridged cysteines lose HG
    bridged = gro(["   57CYS     SG   1   1.000   1.000   1.000",
                   "  146CYS     SG   2   1.200   1.000   1.000",
                   "  146CYS     HG   3   1.300   1.000   1.000"])
    with pytest.raises(ValueError, match="no HG"):
        assert_topology_disulfide_free(clean_top, bridged)

    # 2. force field renamed the residue instead
    cyx_top = tmp_path / "cyx.top"
    cyx_top.write_text("[ atoms ]\n  1  N  57  CYS2  N  1\n")
    with pytest.raises(ValueError, match="re-formed a disulfide"):
        assert_topology_disulfide_free(cyx_top, reduced)


def test_disulfide_prompt_count_covers_all_cysteine_pairs():
    """pdb2gmx asks once per candidate PAIR, not per cysteine.

    SOD1 has 4 cysteines -> up to 6 prompts. Supplying only 4 answers let pdb2gmx take
    its default (form the bond) on the remainder -- the exact invariant being defended.
    """
    from src.fep.pmx_engine import _count_cys_pairs

    pdb = Path(__file__).parent / "fixtures" / "_cys4.pdb"
    pdb.write_text("".join(
        f"ATOM  {i:5d}  SG  CYS A{i:4d}       0.000   0.000   0.000  1.00  0.00           S\n"
        for i in (6, 57, 111, 146)))
    try:
        assert _count_cys_pairs(pdb) == 6      # 4*3/2, not 4
    finally:
        pdb.unlink(missing_ok=True)


def test_mdrun_retries_a_busy_gpu_but_not_a_blown_up_system(monkeypatch, tmp_path):
    """CUDA #46 is transient contention; a LINCS explosion is not.

    Submitting the G93A and I113T arrays in the same instant on 2026-08-09 put their
    first tasks on the same device and lost 4 windows to "CUDA-capable device(s) is/are
    busy". Retrying that is right. Retrying a genuinely broken system would turn a loud
    failure into a slow one, so the match is deliberately narrow.
    """
    from src.fep import pmx_engine

    cfg = {"cluster": {"gpu_retry_attempts": 3, "gpu_retry_delay_s": 0}}
    monkeypatch.setattr(pmx_engine.time, "sleep", lambda _s: None)

    calls = []
    gpu_busy = ("Error while switching to device #0. CUDA error #46 "
                "(cudaErrorDevicesUnavailable): CUDA-capable device(s) is/are busy")

    def flaky(cmd, cwd, stdin=None, dry_run=False, env=None):
        calls.append(cmd)
        if len(calls) < 3:
            raise pmx_engine.ToolError(gpu_busy)
        return "ok"

    monkeypatch.setattr(pmx_engine, "_run", flaky)
    out = pmx_engine._run_mdrun_with_gpu_retry(
        ["gmx", "mdrun"], cwd=tmp_path, cfg=cfg, resume_argv=["-cpi", "prod.cpt"])
    assert out == "ok" and len(calls) == 3

    # a real failure must surface on the first attempt, not be retried
    calls.clear()

    def exploded(cmd, cwd, stdin=None, dry_run=False, env=None):
        calls.append(cmd)
        raise pmx_engine.ToolError("Fatal error:\nLINCS: bond length 9572052.0000")

    monkeypatch.setattr(pmx_engine, "_run", exploded)
    with pytest.raises(pmx_engine.ToolError, match="LINCS"):
        pmx_engine._run_mdrun_with_gpu_retry(
            ["gmx", "mdrun"], cwd=tmp_path, cfg=cfg, resume_argv=[])
    assert len(calls) == 1

    # and it gives up rather than retrying forever
    calls.clear()

    def always_busy(cmd, cwd, stdin=None, dry_run=False, env=None):
        calls.append(cmd)
        raise pmx_engine.ToolError(gpu_busy)

    monkeypatch.setattr(pmx_engine, "_run", always_busy)
    with pytest.raises(pmx_engine.ToolError, match="CUDA error #46"):
        pmx_engine._run_mdrun_with_gpu_retry(
            ["gmx", "mdrun"], cwd=tmp_path, cfg=cfg, resume_argv=[])
    assert len(calls) == 3


def test_mbar_solver_is_run_once_and_refuses_a_non_finite_result(monkeypatch):
    """One solve per (leg, replicate), and an unconverged solve must not become a number.

    pymbar prints "Failed to reach a solution to within tolerance with hybr: trying next
    method" to STDOUT and carries on with a fallback -- it does not raise. Twelve of those
    appeared for a 3-replicate variant, which is 6 legs x 2 because the free energy and
    the overlap matrix were each constructing their own MBAR.
    """
    from src.fep import analyze

    calls = []

    class FakeMBAR:
        def __init__(self, u_kn, N_k):
            calls.append(1)
            # pymbar emits this on STDERR, not stdout -- capturing only stdout reported
            # an empty solver_notes for a run that had five fallbacks.
            print("Failed to reach a solution to within tolerance with hybr: trying next method",
                  file=sys.stderr)

        def compute_free_energy_differences(self):
            n = 3
            return {"Delta_f": np.full((n, n), 2.0), "dDelta_f": np.full((n, n), 0.1)}

        def compute_overlap(self):
            return {"matrix": np.full((3, 3), 1 / 3)}

    import sys
    import types
    fake = types.ModuleType("pymbar")
    fake.MBAR = FakeMBAR
    monkeypatch.setitem(sys.modules, "pymbar", fake)

    out = analyze.solve_leg_mbar(np.zeros((3, 30)), np.array([10, 10, 10]))
    assert len(calls) == 1, "MBAR must be constructed once, not once per quantity"
    assert out["dg_kT"] == 2.0 and out["min_adjacent"] == pytest.approx(1 / 3)
    # the solver's fallback chatter is captured rather than lost to stdout
    assert any("hybr" in n for n in out["solver_notes"])

    class NonFinite(FakeMBAR):
        def compute_free_energy_differences(self):
            return {"Delta_f": np.full((3, 3), np.nan), "dDelta_f": np.full((3, 3), 0.1)}

    fake.MBAR = NonFinite
    with pytest.raises(FloatingPointError, match="non-finite"):
        analyze.solve_leg_mbar(np.zeros((3, 30)), np.array([10, 10, 10]))


def test_resume_guard_refuses_an_unstamped_checkpoint(tmp_path):
    """An UNSTAMPED run directory is not a fresh one -- its provenance is unknown.

    G93A, 2026-08-11: run dirs predating the fingerprint survived a protocol change,
    mdrun resumed their already-complete prod.cpt (both protocols share nsteps), appended
    nothing, and the stale 176-record dhdl.xvg was sliced by the new discard of 500. All
    108 windows were written as (18, 0) carrying the NEW hash. The guard compared only
    when a stamp existed, so it failed open on exactly the case it was written for.
    """
    from src.fep.pmx_engine import ToolError, assert_resumable

    run_dir = tmp_path / "w0_r0"
    run_dir.mkdir()

    # fresh directory, no checkpoint: allowed, and stamped for next time
    assert_resumable(run_dir, "NEW0000")
    assert (run_dir / "protocol.sha").read_text().strip() == "NEW0000"

    # matching stamp + checkpoint: a legitimate resume
    (run_dir / "prod.cpt").write_bytes(b"cpt")
    assert_resumable(run_dir, "NEW0000")

    # stamp disagrees: refuse
    (run_dir / "protocol.sha").write_text("OLD0000\n")
    with pytest.raises(ToolError, match="OLD0000"):
        assert_resumable(run_dir, "NEW0000")

    # NO stamp but a checkpoint exists: the actual G93A case -- must also refuse
    (run_dir / "protocol.sha").unlink()
    with pytest.raises(ToolError, match="pre-fingerprint"):
        assert_resumable(run_dir, "NEW0000")


def test_window_with_no_production_samples_is_refused():
    """(n_states, 0) must raise here, not IndexError inside pymbar much later.

    The real G93A numbers: 176 records in a stale dhdl.xvg against a discard of 500.
    Written silently, it surfaced as `IndexError: index 0 is out of bounds for axis 0
    with size 0` from mbar_solvers, with nothing pointing at the cause.
    """
    from src.fep.pmx_engine import discard_equilibration

    ok = discard_equilibration(np.zeros((18, 3501)), 500, "A4V/folded w0 r0")
    assert ok.shape == (18, 3001)

    with pytest.raises(ValueError, match="176 records"):
        discard_equilibration(np.zeros((18, 176)), 500, "G93A/folded w0 r0")


def test_persistent_gpu_unavailability_asks_sge_to_reschedule(monkeypatch, tmp_path):
    """Waiting cannot fix a device this host will never give us -- only moving can.

    SGE's `gpus` complex is bookkeeping, not enforcement, and the cards are in
    Exclusive_Process mode, so another user's context makes a device unusable while the
    scheduler still counts it free. A4V/G93A/F64A each lost exactly the first task of a
    leg this way, F64A even with 12 retries over 24 minutes on two DIFFERENT nodes.
    Exit 99 is Grid Engine's "reschedule this task".
    """
    from src.fep import pmx_engine

    cfg = {"cluster": {"gpu_retry_attempts": 3, "gpu_retry_delay_s": 0}}
    monkeypatch.setattr(pmx_engine.time, "sleep", lambda _s: None)

    def always_busy(cmd, cwd, stdin=None, dry_run=False, env=None):
        raise pmx_engine.ToolError(
            "Error while switching to device #0. CUDA error #46 "
            "(cudaErrorDevicesUnavailable): CUDA-capable device(s) is/are busy")

    monkeypatch.setattr(pmx_engine, "_run", always_busy)
    with pytest.raises(pmx_engine.GpuUnavailableError, match="place this task elsewhere"):
        pmx_engine._run_mdrun_with_gpu_retry(
            ["gmx", "mdrun"], cwd=tmp_path, cfg=cfg, resume_argv=[])

    # a genuine failure must stay a plain ToolError -- rescheduling it would loop forever
    def exploded(cmd, cwd, stdin=None, dry_run=False, env=None):
        raise pmx_engine.ToolError("Fatal error:\nLINCS: bond length 9572052.0000")

    monkeypatch.setattr(pmx_engine, "_run", exploded)
    with pytest.raises(pmx_engine.ToolError) as caught:
        pmx_engine._run_mdrun_with_gpu_retry(
            ["gmx", "mdrun"], cwd=tmp_path, cfg=cfg, resume_argv=[])
    assert not isinstance(caught.value, pmx_engine.GpuUnavailableError)


def test_solver_notes_drops_the_jax_banner_but_keeps_real_warnings():
    """8800b73 captured stderr so a marginal MBAR solve could not report as clean.

    pymbar's JAX backend then filled solver_notes with its six-line import banner. Nothing
    was lost (the notes are an unbounded set) but a field that is pure boilerplate stops
    being read, which defeats the reason stderr is captured at all.
    """
    from src.fep.analyze import _solver_notes

    banner = (
        "******************************************\n"
        "******* JAX 64-bit mode is now on! *******\n"
        "*   This MAY cause problems with other   *\n"
        "*      uses of JAX in the same code.     *\n"
        "*     JAX is now set to 64-bit mode!     *\n"
        "******************************************"
    )
    assert _solver_notes(banner) == set()

    warning = "Failed to reach a solution to within tolerance with hybr: trying next method"
    assert _solver_notes(banner + "\n" + warning) == {warning}


def test_gate_refuses_to_correlate_variants_from_different_protocols():
    """analyze._check_single_protocol guards WITHIN a variant; nothing guarded ACROSS them.

    A gate could mix a pre-fix ΔΔG with a post-fix one and correlate numbers that are not
    on the same scale. Fails closed -- which protocol is correct is not evaluate_gate's
    decision to make.
    """
    from src.analysis.validate import evaluate_gate

    cfg = {"validation": {"gate_subset": ["A", "B", "C"], "min_pearson": 0.7,
                          "min_gate_points": 2},
           "fep": {"framework": "gromacs_pmx",
                   "convergence": {"max_cycle_closure_kcal": 1.0}}}
    exp = {"A": 1.0, "B": 2.0, "C": 3.0}

    def rec(ddg, protocol):
        return {"ddg": ddg, "cycle_closure_kcal": 0.2, "converged": True,
                "provenance": "gromacs_pmx", "protocol": protocol}

    mixed = {"A": rec(1.1, "aaa"), "B": rec(2.1, "bbb"), "C": rec(3.1, "aaa")}
    verdict = evaluate_gate(mixed, exp, cfg)
    assert verdict["passed"] is False
    assert "multiple protocols" in verdict["reason"]

    same = {v: rec(d + 0.1, "aaa") for v, d in exp.items()}
    passing = evaluate_gate(same, exp, cfg)
    assert passing["passed"] is True
    assert passing["protocol"] == "aaa"


def test_lambda_ladder_defaults_to_uniform_and_validates_an_explicit_vector():
    """fep.lambda_vector moves windows where overlap is thin without resizing the array."""
    import pytest as _pytest

    from src.fep.pmx_engine import lambda_ladder

    uniform = lambda_ladder({"fep": {"lambda_windows": 5}})
    assert uniform == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert lambda_ladder({"fep": {"lambda_windows": 5, "lambda_vector": None}}) == uniform

    dense_high = [0.0, 0.4, 0.7, 0.9, 1.0]
    assert lambda_ladder({"fep": {"lambda_windows": 5, "lambda_vector": dense_high}}) == dense_high

    # Wrong length would desynchronise submit_array.sh's array size and dhdl_to_u_kn.
    with _pytest.raises(ValueError, match="lambda_windows"):
        lambda_ladder({"fep": {"lambda_windows": 18, "lambda_vector": dense_high}})
    with _pytest.raises(ValueError, match="ascend"):
        lambda_ladder({"fep": {"lambda_windows": 5,
                               "lambda_vector": [0.0, 0.4, 0.4, 0.9, 1.0]}})
    with _pytest.raises(ValueError, match="0.0 to 1.0"):
        lambda_ladder({"fep": {"lambda_windows": 5,
                               "lambda_vector": [0.0, 0.4, 0.7, 0.9, 0.95]}})


def test_independent_replicate_systems_changes_the_path_seed_and_protocol():
    """The flag alters every window's starting structure while leaving the .mdp identical.

    Without it in the fingerprint a variant could mix shared-box and independent-box
    replicates and the guard would see one protocol. And with the flag OFF the hash must
    be unchanged, or turning the key on in config invalidates data it does not affect.
    """
    import copy

    import pytest as _pytest
    import yaml

    from src.fep.pmx_engine import (build_system, production_mdp, protocol_extra,
                                    protocol_fingerprint, system_dir)

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    shared = copy.deepcopy(cfg)
    shared["fep"]["independent_replicate_systems"] = False
    per_rep = copy.deepcopy(cfg)
    per_rep["fep"]["independent_replicate_systems"] = True

    # Shared: one directory, replicate ignored. Per-replicate: one directory each.
    assert system_dir(shared, "F64A", "folded", 0) == system_dir(shared, "F64A", "folded", 2)
    assert system_dir(shared, "F64A", "folded", 0).name == "system"
    assert system_dir(per_rep, "F64A", "folded", 0) != system_dir(per_rep, "F64A", "folded", 2)
    assert system_dir(per_rep, "F64A", "folded", 2).name == "system_r2"

    # A per-replicate build with no replicate would hand every replicate one box again.
    with _pytest.raises(ValueError, match="per replicate"):
        system_dir(per_rep, "F64A", "folded", None)
    with _pytest.raises(ValueError, match="per replicate"):
        build_system(per_rep, "F64A", "folded", rep=None, dry_run=True)

    def fp(c):
        return protocol_fingerprint(production_mdp(c, 0, 1, leg="folded"), protocol_extra(c))

    assert protocol_extra(shared) == {}                 # default adds nothing to the hash
    assert fp(shared) == protocol_fingerprint(production_mdp(shared, 0, 1, leg="folded"))
    assert fp(per_rep) != fp(shared)                    # flipping it IS a protocol change


def test_submit_array_task_count_matches_the_config():
    """`#$ -t` is a static SGE directive; the array size lives in config. They drift.

    submit_array.sh re-derives legs*windows*replicates at run time and aborts on a
    mismatch -- but only after the array is queued, so every task exits 2 and the wall
    time is wasted. This catches it at commit time instead. CLAUDE.md rule 5 and the
    lambda_windows 18 -> 20 change on 2026-08-15 are both about exactly this.
    """
    import re

    import yaml

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    fcfg = cfg["fep"]
    expected = len(fcfg["legs"]) * int(fcfg["lambda_windows"]) * int(fcfg["replicates"])

    script = (ROOT / "scripts" / "submit_array.sh").read_text()
    directives = re.findall(r"^#\$ -t 1-(\d+)\s*$", script, flags=re.MULTILINE)
    assert len(directives) == 1, f"expected exactly one '#$ -t' directive, got {directives}"
    assert int(directives[0]) == expected, (
        f"submit_array.sh declares -t 1-{directives[0]} but "
        f"legs*lambda_windows*replicates = {expected}. Update the directive and the "
        "config in the SAME commit."
    )


def test_lambda_vector_in_config_is_the_ladder_the_mdp_gets():
    """The configured vector must survive into fep-lambdas at the .mdp's 4dp precision.

    A vector written to more precision than '%.4f' would round in the .mdp, so the ladder
    actually sampled would differ from the one in config -- and the protocol hash would
    still look fine because it hashes the .mdp.
    """
    import yaml

    from src.fep.pmx_engine import lambda_ladder, production_mdp

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    ladder = lambda_ladder(cfg)
    assert len(ladder) == int(cfg["fep"]["lambda_windows"])

    line = next(l for l in production_mdp(cfg, 0, 1, leg="folded").splitlines()
                if l.startswith("fep-lambdas"))
    emitted = [float(x) for x in line.split("=", 1)[1].split()]
    assert emitted == [round(x, 4) for x in ladder]
    assert emitted == ladder, "config vector rounds when written to the .mdp"


def test_per_leg_sampling_reaches_the_mdp_and_the_discard():
    """fep.ns_per_window may differ per leg, and BOTH consumers must honour it.

    The unfolded tripeptide reproduces to ~0.1 kcal/mol across independent boxes while
    the folded protein spreads to 1.46, so sampling them equally spends about a third of
    the compute where there is no defect. The hazard is that the .mdp and the
    equilibration discard are computed in two places; if they disagree, the discard
    misaligns with what was actually written and the ddG is plausible and wrong.
    """
    import copy
    import yaml
    from src.fep.pmx_engine import leg_value, production_mdp, window_schedule

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    cfg["fep"]["ns_per_window"] = {"folded": 9, "unfolded": 3}
    cfg["fep"]["equilibration_ns"] = {"folded": 2.0, "unfolded": 0.5}

    assert leg_value(cfg, "ns_per_window", "folded") == 9
    assert leg_value(cfg, "ns_per_window", "unfolded") == 3

    f = window_schedule(cfg, "folded")
    u = window_schedule(cfg, "unfolded")
    # 9 ns + 2 ns at 2 fs = 5.5M steps; 3 ns + 0.5 ns = 1.75M
    assert f["nsteps"] == 5_500_000 and u["nsteps"] == 1_750_000
    # the discard must cover the equilibration each leg actually ran
    assert f["discard"] * f["nstdhdl"] >= f["equil_steps"]
    assert u["discard"] * u["nstdhdl"] >= u["equil_steps"]
    assert f["discard"] != u["discard"]

    # and the .mdp the leg gets must carry ITS nsteps, not the other leg's
    fm = production_mdp(cfg, 0, 1, leg="folded")
    um = production_mdp(cfg, 0, 1, leg="unfolded")
    assert "nsteps                   = 5500000" in fm
    assert "nsteps                   = 1750000" in um

    # a scalar config still works and gives both legs the same schedule
    scalar = copy.deepcopy(cfg)
    scalar["fep"]["ns_per_window"] = 3
    scalar["fep"]["equilibration_ns"] = 0.5
    assert window_schedule(scalar, "folded") == window_schedule(scalar, "unfolded")


def test_per_leg_setting_without_a_leg_raises_rather_than_guessing():
    """Silently defaulting would give the two legs one sampling length by accident."""
    import yaml
    from src.fep.pmx_engine import leg_value, window_schedule

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    cfg["fep"]["ns_per_window"] = {"folded": 9, "unfolded": 3}
    with pytest.raises(ValueError, match="per-leg"):
        leg_value(cfg, "ns_per_window", None)
    with pytest.raises(ValueError, match="per-leg"):
        window_schedule(cfg, None)
    with pytest.raises(ValueError, match="no entry for leg"):
        leg_value(cfg, "ns_per_window", "sideways")


def test_per_leg_sampling_makes_the_two_legs_hash_differently():
    """Which is why _check_single_protocol had to become per-leg.

    A leg's windows must agree with each other; the two legs need not agree with one
    another, because each is a separate MBAR estimate and they never share samples.
    """
    import yaml
    from src.fep.pmx_engine import production_mdp, protocol_fingerprint

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    cfg["fep"]["ns_per_window"] = {"folded": 9, "unfolded": 3}
    fp = lambda leg, w=0, seed=1: protocol_fingerprint(
        production_mdp(cfg, w, seed, leg=leg))
    assert fp("folded") != fp("unfolded")          # differ by design
    assert fp("folded") == fp("folded", w=19, seed=999)   # invariant within a leg


def test_run_pmx_window_uses_THIS_leg_schedule_for_both_mdp_and_discard(monkeypatch, tmp_path):
    """The WIRING, not the helpers.

    Two earlier regressions in this project were guards that existed, were unit-tested,
    and were not connected. Testing window_schedule() and production_mdp() in isolation
    passes even if run_pmx_window hardcodes the wrong leg or drops the argument, so this
    drives the real function with GROMACS stubbed out and asserts that a folded window
    gets the FOLDED nsteps and the FOLDED discard.
    """
    import yaml
    from src.fep import pmx_engine as pe

    cfg = yaml.safe_load(open(ROOT / "config" / "pipeline.yaml"))
    cfg["fep"]["ns_per_window"] = {"folded": 9, "unfolded": 3}
    cfg["fep"]["equilibration_ns"] = {"folded": 2.0, "unfolded": 0.5}

    written: dict[str, str] = {}
    monkeypatch.setattr(pe, "build_system", lambda *a, **k: tmp_path / "system")
    monkeypatch.setattr(pe, "gmx_command", lambda cfg: "gmx")
    monkeypatch.setattr(pe, "_run", lambda *a, **k: "")
    monkeypatch.setattr(pe, "_run_mdrun_with_gpu_retry", lambda *a, **k: "")
    monkeypatch.setattr(pe, "assert_resumable", lambda *a, **k: None)
    monkeypatch.setattr(pe, "ROOT", tmp_path)

    def fake_write(path, body, dry_run=False):
        written[Path(path).name] = body
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(body)
        return Path(path)
    monkeypatch.setattr(pe, "_write_mdp", fake_write)

    # a dhdl with one record per nstdhdl over the FOLDED nsteps, plus t=0
    f_sched = pe.window_schedule(cfg, "folded")
    n_records = f_sched["nsteps"] // f_sched["nstdhdl"] + 1
    monkeypatch.setattr(pe, "dhdl_to_u_kn",
                        lambda *a, **k: np.zeros((cfg["fep"]["lambda_windows"], n_records)))

    (tmp_path / "results" / "fep" / "A4V" / "folded" / "w0_r0").mkdir(parents=True)
    (tmp_path / "results" / "fep" / "A4V" / "folded" / "w0_r0" / "em.gro").write_text("x")

    out = pe.run_pmx_window(cfg, "A4V", "folded", 0, 0)

    assert "nsteps                   = 5500000" in written["prod.mdp"], \
        "folded window did not get the folded sampling length"
    # discard must match the equilibration THIS leg actually ran, not the other leg's
    assert out["u_kn_window"].shape[1] == n_records - f_sched["discard"]
    u_sched = pe.window_schedule(cfg, "unfolded")
    assert f_sched["discard"] != u_sched["discard"]
