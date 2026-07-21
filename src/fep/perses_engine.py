"""Real OpenMM+Perses alchemical FEP window (GPU-only; runs on the BU SCC).

This is the science kernel behind :func:`src.fep.window._perses_window`. One call runs ONE
(variant, leg, window, replicate): it builds the WT->mut hybrid topology with Perses,
places it at this window's lambda on the validated default ``LambdaProtocol``, runs Langevin
dynamics, and scores every collected frame at ALL lambda states -> ``u_kn_window`` for MBAR
(the schema :mod:`src.fep.analyze` already consumes).

Design decisions (see memory ``fep-window-design``):
  * folded leg  -> WT apo/reduced protein (from ``prep.build.prepare_wt_apo``), Perses solvates;
  * unfolded leg-> capped tripeptide with native neighbours (``prep.tripeptide``), Perses solvates;
  * MONOMER only for now -- the dimer double-mutation case raises (deferred to M4);
  * one window == one GPU job, so we sample a single lambda per process (no built-in repex)
    and cross-evaluate energies at the other states ourselves.

Everything numeric is read from ``config/pipeline.yaml`` (CLAUDE.md #4). This module cannot
run without a GPU + the SCC ``sod1-fep`` env (openmmtools/perses), so imports of those are
lazy and the pure/guard logic stays importable for tests.

NOTE (perses API): kwarg names on ``PointMutationExecutor`` and the alchemical-state helpers
have drifted across perses releases. This targets perses 0.10.x; if the SCC smoke-test raises
a TypeError on a kwarg, that constructor call is the place to adjust -- the surrounding
sampling/MBAR logic is engine-stable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from src.prep.build import load_variant_record

ROOT = Path(__file__).resolve().parents[2]

# config value -> OpenMM forcefield xml. Explicit map so an unknown value fails loudly
# (CLAUDE.md #4) instead of silently defaulting to the wrong force field.
_FF_PROTEIN = {"amber14sb": "amber14/protein.ff14SB.xml"}
_FF_WATER = {
    "tip3p": "amber14/tip3p.xml",
    "tip3pfb": "amber14/tip3pfb.xml",
    "tip4pew": "amber14/tip4pew.xml",
}
_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def _forcefield_files(cfg: dict) -> list[str]:
    scfg = cfg["structure"]
    try:
        return [_FF_PROTEIN[scfg["ff_protein"]], _FF_WATER[scfg["water_model"]]]
    except KeyError as e:  # surface the gap rather than inventing a default
        raise KeyError(
            f"No force-field xml mapped for {e!s}; add it to perses_engine._FF_* "
            "(structure.ff_protein / structure.water_model)."
        ) from e


def _lambda_value(window: int, n_states: int) -> float:
    """Evenly spaced global lambda in [0, 1] for ``n_states`` windows."""
    return window / (n_states - 1)


def _window_geometry(cfg: dict, smoke: bool) -> dict:
    """Resolve step counts / frame cadence from config (or a tiny smoke-test override)."""
    fcfg = cfg["fep"]
    dt_fs = float(fcfg["timestep_fs"])
    steps_per_ns = int(round(1_000_000.0 / dt_fs))  # 1 ns = 1e6 fs
    if smoke:  # cheap first-light run to shake out the GPU/Perses path
        return {"dt_fs": dt_fs, "equil_steps": 100, "frames": 4, "steps_per_frame": 50}
    frames = int(fcfg["frames_per_window"])
    prod_steps = int(round(float(fcfg["ns_per_window"]) * steps_per_ns))
    return {
        "dt_fs": dt_fs,
        "equil_steps": int(round(float(fcfg["equilibration_ns"]) * steps_per_ns)),
        "frames": frames,
        "steps_per_frame": max(1, prod_steps // frames),
    }


def _input_structure(cfg: dict, variant: str, leg: str) -> Path:
    """Build (once, memoised) and return the WT input structure Perses consumes for ``leg``.

    Cached under ``results/fep/<variant>/inputs/`` so the 180 array tasks of a variant reuse
    one build. Written atomically (temp + rename) to tolerate concurrent array tasks.
    """
    from src.prep.build import load_config as _lc  # noqa: F401  (keep import graph explicit)
    from src.prep.build import prepare_wt_apo
    from src.prep.tripeptide import build_tripeptide

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    if leg == "folded" and record["oligomer"] == "dimer":
        raise NotImplementedError(
            f"{variant} is a dimer; the folded leg would carry the mutation in both subunits "
            "(double alchemical transformation). Deferred to M4 -- monomer only for the gate."
        )

    in_dir = ROOT / "results" / "fep" / variant / "inputs"
    in_dir.mkdir(parents=True, exist_ok=True)

    if leg == "folded":
        target = in_dir / "wt_apo.pdb"
        if not target.exists():
            tmp = target.with_suffix(".pdb.tmp")
            prepare_wt_apo(cfg, variant, tmp)
            tmp.replace(target)
        return target
    if leg == "unfolded":
        wt = _input_structure(cfg, variant, "folded")  # tripeptide neighbours come from WT
        target = in_dir / "tripeptide.pdb"
        if not target.exists():
            tmp = target.with_suffix(".pdb.tmp")
            build_tripeptide(cfg, variant, wt, tmp)
            tmp.replace(target)
        return target
    raise ValueError(f"unknown leg {leg!r} (expected folded|unfolded)")


def _build_htf(cfg: dict, variant: str, leg: str):
    """Build the WT->mut hybrid topology factory for this (variant, leg) via Perses."""
    from openmm import app, unit
    from perses.app.relative_point_mutation_setup import PointMutationExecutor

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    scfg, fcfg = cfg["structure"], cfg["fep"]
    pdb = _input_structure(cfg, variant, leg)

    constraints = {"HBonds": app.HBonds, "AllBonds": app.AllBonds, "None": None}[fcfg["constraints"]]
    ff_kwargs = {
        "removeCMMotion": False,
        "constraints": constraints,
        "rigidWater": True,
        "hydrogenMass": (float(fcfg["hydrogen_mass_amu"]) * unit.amu)
        if fcfg.get("hydrogen_mass_amu") else None,
    }
    periodic_ff_kwargs = {
        "nonbondedMethod": {"PME": app.PME}[fcfg["nonbonded_method"]],
        "nonbondedCutoff": float(fcfg["nonbonded_cutoff_nm"]) * unit.nanometer,
    }

    pm = PointMutationExecutor(
        protein_filename=str(pdb),
        mutation_chain_id="A",                       # monomer ref chain (guarded above)
        mutation_residue_id=str(record["mature_pos"]),
        proposed_residue=_THREE[record["mut_aa"]],
        forcefield_files=_forcefield_files(cfg),
        water_model=scfg["water_model"],
        ionic_strength=float(scfg["ion_conc_M"]) * unit.molar,
        solvent_padding=float(scfg["solvent_padding_nm"]) * unit.nanometer,
        temperature=float(fcfg["temperature_K"]) * unit.kelvin,
        forcefield_kwargs=ff_kwargs,
        periodic_forcefield_kwargs=periodic_ff_kwargs,
        conduct_endstate_validation=False,
    )
    return pm.get_apo_htf()  # apo phase = solvated (protein or tripeptide), no ligand


def _compound_states(cfg: dict, hybrid_system):
    """One CompoundThermodynamicState per lambda window, along the default LambdaProtocol."""
    from openmm import unit
    from openmmtools.states import CompoundThermodynamicState, ThermodynamicState
    from perses.annihilation.lambda_protocol import LambdaProtocol, RelativeAlchemicalState

    fcfg = cfg["fep"]
    n_states = int(fcfg["lambda_windows"])
    T = float(fcfg["temperature_K"]) * unit.kelvin
    P = (float(fcfg["pressure_atm"]) * unit.atmosphere) if fcfg["ensemble"] == "NPT" else None
    protocol = LambdaProtocol(functions=fcfg["lambda_protocol"])

    states = []
    for w in range(n_states):
        alch = RelativeAlchemicalState.from_system(hybrid_system)
        thermo = ThermodynamicState(hybrid_system, temperature=T, pressure=P)
        cs = CompoundThermodynamicState(thermo, composable_states=[alch])
        cs.set_alchemical_parameters(_lambda_value(w, n_states), protocol)
        states.append(cs)
    return states


def run_perses_window(cfg: dict, variant: str, leg: str, window: int, rep: int,
                      smoke: bool = False) -> dict:
    """Sample one lambda window on GPU and return reduced potentials at every state.

    Returns the ``src.fep.window`` NPZ dict: ``lambda_index``, ``n_states``, ``u_kn_window``
    of shape ``(n_states, n_frames)``.
    """
    from openmm import LangevinMiddleIntegrator, LocalEnergyMinimizer, Platform, unit
    from openmmtools.states import SamplerState

    fcfg = cfg["fep"]
    n_states = int(fcfg["lambda_windows"])
    T = float(fcfg["temperature_K"]) * unit.kelvin
    friction = float(fcfg["collision_rate_per_ps"]) / unit.picosecond
    geom = _window_geometry(cfg, smoke)
    dt = geom["dt_fs"] * unit.femtosecond

    htf = _build_htf(cfg, variant, leg)
    states = _compound_states(cfg, htf.hybrid_system)
    sampler_state = SamplerState(
        positions=htf.hybrid_positions,
        box_vectors=htf.hybrid_system.getDefaultPeriodicBoxVectors(),
    )

    # dynamics context pinned at THIS window's state, on the scheduler-assigned GPU
    # (CUDA_VISIBLE_DEVICES already masks the device -- never set DeviceIndex, CLAUDE.md).
    platform = Platform.getPlatformByName("CUDA")
    integrator = LangevinMiddleIntegrator(T, friction, dt)
    context = states[window].create_context(integrator, platform)
    sampler_state.apply_to_context(context)
    if fcfg.get("minimize", True):
        LocalEnergyMinimizer.minimize(context)
    context.setVelocitiesToTemperature(T)

    ckpt = ROOT / "results" / "fep" / variant / leg / f"w{window}_r{rep}.ckpt.npz"
    u_kn = np.zeros((n_states, geom["frames"]))
    start_frame = 0
    if ckpt.exists() and not smoke:                       # resume a killed window
        start_frame, u_kn, sampler_state = _load_checkpoint(ckpt, u_kn, sampler_state, context)

    if start_frame == 0:
        integrator.step(geom["equil_steps"])

    ckpt_every = _frames_between_checkpoints(cfg, geom)
    for n in range(start_frame, geom["frames"]):
        integrator.step(geom["steps_per_frame"])
        sampler_state.update_from_context(context)
        for l in range(n_states):                          # cross-evaluate at all states
            states[l].apply_to_context(context)
            u_kn[l, n] = states[l].reduced_potential(context)
        states[window].apply_to_context(context)           # restore dynamics state
        if not smoke and (n + 1) % ckpt_every == 0:
            _save_checkpoint(ckpt, n + 1, u_kn, sampler_state)

    ckpt.unlink(missing_ok=True)
    return {"lambda_index": window, "n_states": n_states, "u_kn_window": u_kn}


# --------------------------- checkpoint/resume (windows die) ---------------------------

def _frames_between_checkpoints(cfg: dict, geom: dict) -> int:
    steps = max(1, int(cfg["fep"].get("checkpoint_interval_steps", 50000)))
    return max(1, steps // geom["steps_per_frame"])


def _save_checkpoint(path: Path, next_frame: int, u_kn: np.ndarray, sampler_state) -> None:
    from openmm import unit
    pos = sampler_state.positions.value_in_unit(unit.nanometer)
    box = sampler_state.box_vectors.value_in_unit(unit.nanometer)
    tmp = path.with_suffix(".npz.tmp")
    np.savez(tmp, next_frame=next_frame, u_kn=u_kn, positions=pos, box=box)
    tmp.replace(path)


def _load_checkpoint(path: Path, u_kn: np.ndarray, sampler_state, context):
    from openmm import unit
    d = np.load(path)
    u_kn[:] = d["u_kn"]
    sampler_state.positions = d["positions"] * unit.nanometer
    sampler_state.box_vectors = d["box"] * unit.nanometer
    sampler_state.apply_to_context(context)
    context.setVelocitiesToTemperature(context.getIntegrator().getTemperature())
    return int(d["next_frame"]), u_kn, sampler_state
