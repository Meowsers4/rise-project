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
    """Same engine, different .mdp settings, is still unmixable.

    Provenance catches a mock or another engine. It cannot see sc-coul flipped or nstdhdl
    changed halfway through a variant -- which is exactly what happens if the SCC is
    `git pull`ed while an array is in flight.
    """
    from src.fep.analyze import _check_single_protocol

    _check_single_protocol("A4V", {"abc123"})            # one protocol: fine
    with pytest.raises(ValueError, match="different protocols"):
        _check_single_protocol("A4V", {"abc123", "def456"})
    with pytest.raises(ValueError, match="different protocols"):
        _check_single_protocol("A4V", {"abc123", "unlabelled"})


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
        return protocol_fingerprint(production_mdp(c, window, seed))

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
