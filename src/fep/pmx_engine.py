"""Real GROMACS + pmx alchemical FEP window (GPU; runs on the BU SCC).

The science kernel behind :func:`src.fep.window._pmx_window`. One call runs ONE
(variant, leg, window, replicate): it builds the WT->mut hybrid topology with pmx, runs
that window's lambda state under GROMACS, and converts the resulting ``dhdl.xvg`` into
the ``u_kn_window`` matrix :mod:`src.fep.analyze` already consumes -- the same NPZ schema
the perses engine emitted, so analysis, provenance and the gate are unchanged.

Pipeline, per (variant, leg) -- built ONCE and reused by every window/replicate:

1. ``pmx mutate``   -- WT structure -> hybrid structure carrying both end states
2. ``pdb2gmx``      -- hybrid structure -> .gro + .top under the pmx mutation force field
3. ``pmx gentop``   -- rewrite the topology with B-state (mutant) parameters
4. ``editconf`` / ``solvate`` / ``genion`` -- box, water, neutralising ions
5. energy minimisation

then per window: write the ``.mdp`` with this window's ``init-lambda-state``, ``grompp``,
``mdrun``, and parse ``dhdl.xvg``.

Two project invariants need explicit defence here, because GROMACS would otherwise undo
Stage 1's work:

* **Reduced disulfide.** ``prep.build`` breaks Cys57-Cys146 on an OpenMM topology, but
  ``pdb2gmx`` re-detects disulfides from SG-SG distance and would rebuild it. We feed it
  ``-ss`` answers that keep every cysteine a free thiol (``fep.keep_disulfide_reduced``).
* **One Hamiltonian per leg.** All windows of a leg must sample the same system, so the
  solvated/minimised system is built once, cached, and reused -- never rebuilt per job.

Everything numeric comes from ``config/pipeline.yaml`` (CLAUDE.md #4). NOTE: the SCC's
GROMACS is an OpenMPI build, so ``mdrun`` runs under ``cluster.mdrun_launcher`` and
``-ntmpi`` is unavailable; threads come from ``cluster.mdrun_ntomp``.
"""

from __future__ import annotations

import fcntl
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from src.prep.build import load_variant_record
from src.seeds import stable_seed

ROOT = Path(__file__).resolve().parents[2]

# Written by _build_system() when a leg's system is complete; its absence means a previous
# build died partway and the directory must not be trusted.
_DONE_MARKER = "SYSTEM_READY"

# Candidate force-field directory names inside ``pmx/data``, most specific first. See
# :func:`pmx_ff_dir` -- pmx renamed this between the 2.0 snapshots and current develop.
_PMX_FF_DIRNAMES = ("mutff45", "mutff")


def _log(msg: str) -> None:
    """Progress to stderr -- survives on the SGE task log even if the job is killed."""
    print(f"[pmx_engine] {msg}", file=sys.stderr, flush=True)


class ToolError(RuntimeError):
    """A GROMACS/pmx command failed; carries the tail of its output."""


def _run(cmd: list[str], cwd: Path, stdin: str | None = None, dry_run: bool = False,
         env: dict | None = None) -> str:
    """Run one external command, logging it, and fail loudly with its output.

    Args:
        cmd: argv list.
        cwd: working directory (every stage runs inside the system dir).
        stdin: text piped to the process -- pdb2gmx's interactive prompts need this.
        dry_run: print the command and return without executing.
        env: environment overrides merged onto ``os.environ``.

    Returns:
        Combined stdout+stderr.

    Raises:
        ToolError: If the command exits non-zero.
    """
    printable = " ".join(cmd) + (f"   <<< {stdin!r}" if stdin else "")
    _log(f"$ {printable}")
    if dry_run:
        return ""
    t0 = time.time()
    proc = subprocess.run(
        cmd, cwd=str(cwd), input=stdin, capture_output=True, text=True,
        env={**os.environ, **(env or {})},
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        hint = ""
        # exit 127 / a loader error means the binary exists but its runtime libraries
        # (MKL, CUDA) are absent -- i.e. the modules never loaded, not a GROMACS problem.
        if proc.returncode == 127 or "error while loading shared libraries" in out:
            hint = (
                "\n\nHINT: the binary was found but could not load its libraries -- the "
                "environment modules did not load.\n"
                "      Run `source scripts/scc_env.sh` and check it does not warn about "
                "'module: command not found'.\n"
                "      In a non-interactive shell (qrsh ... bash -lc '...') Lmod may be "
                "uninitialised; scc_env.sh bootstraps it, so source it FIRST."
            )
        raise ToolError(
            f"command failed (exit {proc.returncode}) in {cwd}:\n  {printable}\n"
            f"--- last 40 lines ---\n" + "\n".join(out.splitlines()[-40:]) + hint
        )
    _log(f"  ok ({time.time() - t0:.1f}s)")
    return out


# --------------------------------------------------------------------------- #
# Tool resolution                                                              #
# --------------------------------------------------------------------------- #
def gmx_command(cfg: dict) -> str:
    """The GROMACS binary name (``cluster.gmx_binary``; auto-detect if blank)."""
    hint = (
        "Load the GROMACS modules first -- from the repo root:\n"
        "    source scripts/scc_env.sh\n"
        "(that loads cluster.gromacs_prereq_modules then cluster.gromacs_module and "
        "activates the conda env, all from config)."
    )
    name = (cfg["cluster"].get("gmx_binary") or "").strip()
    if name:
        # Resolve even when pinned: otherwise the failure surfaces later as a bare
        # FileNotFoundError from subprocess, which says nothing about modules.
        if shutil.which(name):
            return name
        raise ToolError(f"cluster.gmx_binary={name!r} is not on PATH.\n{hint}")
    for cand in ("gmx", "gmx_mpi"):
        if shutil.which(cand):
            return cand
    raise ToolError(f"no GROMACS binary (gmx / gmx_mpi) on PATH.\n{hint}")


def mdrun_argv(cfg: dict, gpu: bool = True) -> list[str]:
    """argv prefix for mdrun: optional launcher, binary, thread count, GPU offload.

    ``cluster.mdrun_launcher`` is empty for the SCC's serial CUDA build (no MPI, no
    thread-MPI) and would be ``mpirun -np 1`` for the MPI build. ``gpu=False`` drops the
    offload flags for CPU-only steps such as minimisation.
    """
    ccfg = cfg["cluster"]
    argv = [*(ccfg.get("mdrun_launcher") or "").split(), gmx_command(cfg), "mdrun",
            "-ntomp", str(ccfg["mdrun_ntomp"])]
    if gpu:
        argv += (ccfg.get("mdrun_gpu_flags") or "").split()
    return argv


def pmx_argv(cfg: dict, subcommand: str) -> list[str]:
    """argv for a pmx subcommand, preferring the console script over ``-m``.

    pmx's py3 branch exposes a ``pmx`` entry point; if it is absent (e.g. the package was
    put on PYTHONPATH from a clone rather than installed) fall back to the module form.
    """
    if shutil.which("pmx"):
        return ["pmx", subcommand]
    return [sys.executable, "-m", f"pmx.scripts.{subcommand}"]


def pmx_ff_dir(cfg: dict) -> Path:
    """Directory holding pmx's ``<name>.ff`` mutation force fields.

    ``fep.pmx_ff_dir: auto`` derives it from the installed pmx package rather than
    hardcoding a site path, so the config survives a different pmx location.

    The layout is VERSION DEPENDENT and both names are live in the wild:

    * pmx 2.0-era snapshots ship ``data/mutff45`` (a complete GROMACS ff set) alongside a
      legacy flat GROMACS-4.5-style ``data/mutff``; there, mutff45 is the one to use.
    * current develop ships only ``data/mutff``, which is now the maintained set and does
      contain ``amber99sb-star-ildn-mut.ff``.

    So probe both and return whichever actually holds ``fep.pmx_forcefield``. Falling back
    to the first directory that merely exists keeps :func:`verify_tools`' "no such force
    field, available: ..." message useful instead of raising a bare path error here.
    """
    configured = str(cfg["fep"].get("pmx_ff_dir", "auto")).strip()
    if configured and configured != "auto":
        return Path(configured)
    import pmx
    data = Path(pmx.__file__).parent / "data"
    candidates = [data / name for name in _PMX_FF_DIRNAMES]
    ff = cfg["fep"]["pmx_forcefield"]
    for candidate in candidates:
        if (candidate / f"{ff}.ff").is_dir():
            return candidate
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return candidates[0]


def gmx_env(cfg: dict) -> dict:
    """Environment additions for GROMACS: GMXLIB so pdb2gmx sees pmx's force fields."""
    return {"GMXLIB": str(pmx_ff_dir(cfg))}


def verify_tools(cfg: dict) -> dict:
    """Check GROMACS and pmx are callable and the mutation force field exists.

    Returns a dict of what was found. Call this BEFORE submitting an array -- a missing
    module surfaces here in seconds instead of 864 times in the scheduler.
    """
    found = {}
    gmx = gmx_command(cfg)
    found["gmx"] = shutil.which(gmx) or gmx
    version = _run([gmx, "--version"], cwd=ROOT)
    for line in version.splitlines():
        if "GROMACS version" in line or "GPU support" in line:
            found[line.split(":")[0].strip()] = line.split(":", 1)[-1].strip()

    import pmx  # noqa: F401  -- import error here is the answer
    ff = cfg["fep"]["pmx_forcefield"]
    ff_dir = pmx_ff_dir(cfg)
    if not ff_dir.is_dir():
        raise ToolError(f"fep.pmx_ff_dir resolves to {ff_dir}, which does not exist.")
    available = sorted(p.name[:-3] for p in ff_dir.glob("*.ff"))
    if ff not in available:
        raise ToolError(
            f"pmx has no mutation force field {ff!r} in {ff_dir}. Available: {available}. "
            "Set fep.pmx_forcefield to one of those."
        )
    # pdb2gmx must ALSO resolve this force field, which it only does via GMXLIB.
    found["pmx_ff_dir"] = str(ff_dir)
    found["forcefield"] = ff
    found["available_ff"] = ", ".join(available)
    found["GMXLIB"] = gmx_env(cfg)["GMXLIB"]
    return found


# --------------------------------------------------------------------------- #
# System construction (once per variant+leg)                                   #
# --------------------------------------------------------------------------- #
def _input_structure(cfg: dict, variant: str, leg: str) -> Path:
    """Build (once, memoised) the WT input structure for ``leg``.

    Reuses Stage 1 unchanged: ``prepare_wt_apo`` for the folded leg, ``build_tripeptide``
    for the unfolded reference. Both emit PDB, which pdb2gmx consumes directly.
    """
    from src.prep.build import prepare_wt_apo
    from src.prep.tripeptide import build_tripeptide

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    if leg == "folded" and record["oligomer"] == "dimer":
        raise NotImplementedError(
            f"{variant} is a dimer; the folded leg would carry the mutation in both "
            "subunits (double alchemical transformation). Monomer only for the gate."
        )

    in_dir = ROOT / "results" / "fep" / variant / "inputs"
    in_dir.mkdir(parents=True, exist_ok=True)
    lock_path = in_dir / ".inputs.lock"

    def materialize_once(target: Path, builder) -> Path:
        """Build one shared input atomically when folded/unfolded tasks overlap."""
        with open(lock_path, "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if not target.exists():
                tmp = target.with_suffix(target.suffix + ".tmp")
                builder(tmp)
                tmp.replace(target)
        return target

    if leg == "folded":
        target = in_dir / "wt_apo.pdb"
        return materialize_once(target, lambda tmp: prepare_wt_apo(cfg, variant, tmp))
    if leg == "unfolded":
        wt = _input_structure(cfg, variant, "folded")
        target = in_dir / "tripeptide.pdb"
        return materialize_once(target, lambda tmp: build_tripeptide(cfg, variant, wt, tmp))
    raise ValueError(f"unknown leg {leg!r} (expected folded|unfolded)")


def _pdb2gmx_stdin(cfg: dict, n_cysteines: int) -> str:
    """Answers for pdb2gmx's interactive prompts.

    With ``-ss`` interactive, pdb2gmx asks once per cysteine PAIR whether to form a
    disulfide. v1 is the REDUCED form, so every answer is "n" -- this is the GROMACS-side
    guard for the invariant ``prep.build.assert_disulfide_reduced`` enforces upstream.
    """
    if not cfg["fep"].get("keep_disulfide_reduced", True):
        return ""
    return "\n".join(["n"] * max(1, n_cysteines)) + "\n"


def _count_cys_pairs(pdb: Path) -> int:
    """Upper bound on pdb2gmx disulfide prompts: number of SG atoms."""
    return sum(1 for line in pdb.read_text().splitlines()
               if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "SG")


_CAPS = frozenset({"ACE", "NME", "NAC", "NH2"})
_NON_PROTEIN = frozenset({"HOH", "WAT", "SOL", "NA", "CL", "CU", "ZN"})


def _protein_residues(pdb: Path) -> list[tuple[int, str]]:
    """[(resSeq, resname)] for protein residues in file order, caps and solvent dropped."""
    seen: dict[int, str] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[17:20].strip()
        if name in _CAPS or name in _NON_PROTEIN:
            continue
        seen.setdefault(int(line[22:26]), name)
    return list(seen.items())


# pdb2gmx writes FORCE-FIELD residue names, not PDB ones: under amber99sb-star-ildn-mut
# a histidine is HID/HIE/HIP depending on protonation, never HIS. Comparing those against
# the panel's one-letter code would reject every HIS variant in the panel (H43R, H46R,
# H71Q, H110Y, H120L) with a spurious "refusing to mutate the wrong residue".
_RESNAME_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
    "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS",
}


def _canonical_resname(name: str | None) -> str | None:
    """Map a force-field protonation/bonding variant back to its PDB residue name."""
    return _RESNAME_ALIASES.get(name, name) if name else name


def mutation_resid(cfg: dict, variant: str, pdb: Path, leg: str) -> int:
    """Residue id pmx must mutate, READ FROM the structure pmx will actually parse.

    The two legs number residues differently and neither can be assumed:

    * **folded** -- the protein keeps mature numbering, so the panel's ``mature_pos`` is
      the id; we still verify the residue there is the expected wild type.
    * **unfolded** -- the capped tripeptide is renumbered by PDBFixer/pdb2gmx (ACE-X-Y-Z-NME
      becomes 1..5), so the mutation site is the MIDDLE protein residue, not ``mature_pos``.
      Using mature_pos here silently mutates a neighbour and yields a plausible, wrong ΔΔG.

    Raises:
        ValueError: If the residue found is not the panel's wild type -- the guard that
            makes a numbering mistake loud instead of silent.
    """
    from src.prep.build import _THREE

    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    want = _THREE[record["wt_aa"]]
    residues = _protein_residues(pdb)
    if not residues:
        raise ValueError(f"{pdb}: no protein residues parsed")

    if leg == "folded":
        resid = int(record["mature_pos"])
        got = dict(residues).get(resid)
    elif leg == "unfolded":
        if len(residues) != 3:
            raise ValueError(
                f"{pdb}: unfolded reference has {len(residues)} protein residues, "
                "expected the 3 of a capped tripeptide."
            )
        resid, got = residues[1]            # the middle residue is the mutation site
    else:
        raise ValueError(f"unknown leg {leg!r}")

    if _canonical_resname(got) != want:
        raise ValueError(
            f"{pdb}: residue {resid} is {got!r}, panel says {want!r} "
            f"({record['variant']}). Refusing to mutate the wrong residue."
        )
    return resid


def write_mutation_script(cfg: dict, variant: str, pdb: Path, leg: str, path: Path,
                          dry_run: bool = False) -> str:
    """Write pmx's ``-script`` file: one ``"<resid> <target>"`` pair.

    ``pmx mutate -h``: *"The script file simply has to consist of 'resi_number
    target_residue.' pairs"*, in an extended one-letter code. One line, since v1 is a
    single point mutation. The resid comes from :func:`mutation_resid`, i.e. from the
    actual file, never assumed.
    """
    record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
    resid = int(record["mature_pos"]) if dry_run else mutation_resid(cfg, variant, pdb, leg)
    line = f"{resid} {record['mut_aa']}\n"
    _log(f"mutation script ({leg}): {line.strip()!r} -- {record['wt_aa']}"
         f"{record['mature_pos']}{record['mut_aa']} is residue {resid} in {pdb.name}")
    if not dry_run:
        path.write_text(line)
    return line


def build_system(cfg: dict, variant: str, leg: str, dry_run: bool = False) -> Path:
    """Build (once) the solvated, minimised hybrid system for ``variant``/``leg``.

    Cached under ``results/fep/<variant>/<leg>/system/`` and marked complete only when
    every stage succeeded, so a half-built directory is never silently reused. Every
    window and replicate of this leg shares it -- MBAR requires one Hamiltonian.

    All array tasks of a (variant, leg) start at once, so the build is serialised with
    an ``flock``: the first task builds the system, the rest block on the lock and then
    reuse the cached one. Without the lock two tasks can rmtree each other's half-built
    system mid-build.
    """
    fcfg, scfg = cfg["fep"], cfg["structure"]
    sys_dir = ROOT / "results" / "fep" / variant / leg / "system"

    lock_path = sys_dir.parent / ".system.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    # Production grompp also needs GMXLIB when this process reuses an existing system;
    # do this before the marker fast path, not only while constructing the cache.
    os.environ.update(gmx_env(cfg))
    with open(lock_path, "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (sys_dir / _DONE_MARKER).exists():
            _log(f"reusing cached system {sys_dir}")
            return sys_dir

        if sys_dir.exists():                   # previous build died partway
            _log(f"discarding incomplete system dir {sys_dir}")
            shutil.rmtree(sys_dir)
        sys_dir.mkdir(parents=True, exist_ok=True)

        record = load_variant_record(ROOT / cfg["panel"]["csv"], variant)
        gmx = gmx_command(cfg)
        ff = fcfg["pmx_forcefield"]
        water = scfg["water_model"]
        wt_pdb = _input_structure(cfg, variant, leg)
        shutil.copy(wt_pdb, sys_dir / "wt.pdb")

        # Deterministic ion placement: same system for every window/replicate.
        seed = stable_seed(variant, leg, "system")

        n_sg = _count_cys_pairs(sys_dir / "wt.pdb") if not dry_run else 4
        ss_answers = _pdb2gmx_stdin(cfg, n_sg)
        # Hydrogen-mass repartitioning is a pdb2gmx flag in GROMACS, not an .mdp setting.
        # Off by default: a >2 fs timestep needs it, but combining it with pmx's dummy
        # masses is untested here (fep.heavyh).
        heavyh = ["-heavyh"] if fcfg.get("heavyh") else []

        # 1. pdb2gmx FIRST, on the wild type. `pmx mutate -h`: "The best way to use this
        #    script is to take a pdb/gro file that has been written with pdb2gmx with all
        #    hydrogen atoms present." Feeding it a raw PDBFixer file makes pmx's reader
        #    parse zero residues and die in make_chains(). This pass only normalises
        #    naming and hydrogens to the pmx force field; its topology is discarded.
        _run([gmx, "pdb2gmx", "-f", "wt.pdb", "-o", "wt_gmx.pdb", "-p", "wt_discard.top",
              "-ff", ff, "-water", water, "-ignh", "-ss", *heavyh],
             cwd=sys_dir, stdin=ss_answers, dry_run=dry_run)

        # 2. pmx mutate on the pdb2gmx output. The resid is read from THAT file: pdb2gmx
        #    and PDBFixer renumber, so for the tripeptide the site is the middle residue,
        #    not mature_pos (mutation_resid also verifies the wild-type identity there).
        write_mutation_script(cfg, variant, sys_dir / "wt_gmx.pdb", leg,
                              sys_dir / "mutation.txt", dry_run=dry_run)
        _run([*pmx_argv(cfg, "mutate"),
              "-f", "wt_gmx.pdb", "-o", "hybrid.pdb", "-ff", ff,
              "--script", "mutation.txt"],
             cwd=sys_dir, dry_run=dry_run)

        # 3. pdb2gmx again, now on the hybrid -- this is the topology we keep.
        #    No -ignh here: the hybrid residue's hydrogens were built by pmx mutate and
        #    pdb2gmx has no hdb entry for hybrid residues, so rebuilding them fails.
        _run([gmx, "pdb2gmx", "-f", "hybrid.pdb", "-o", "conf.gro", "-p", "topol.top",
              "-ff", ff, "-water", water, "-ss", *heavyh],
             cwd=sys_dir, stdin=ss_answers, dry_run=dry_run)

        # 4. pmx gentop -- write the B-state (mutant) parameters into the topology.
        _run([*pmx_argv(cfg, "gentop"),
              "-p", "topol.top", "-o", "hybrid.top", "-ff", ff],
             cwd=sys_dir, dry_run=dry_run)

        # 4. box, solvent, neutralising ions at the configured ionic strength.
        _run([gmx, "editconf", "-f", "conf.gro", "-o", "box.gro", "-bt", "cubic",
              "-d", str(scfg["solvent_padding_nm"]), "-c"], cwd=sys_dir, dry_run=dry_run)
        _run([gmx, "solvate", "-cp", "box.gro", "-cs", "spc216.gro",
              "-o", "solv.gro", "-p", "hybrid.top"], cwd=sys_dir, dry_run=dry_run)
        _write_mdp(sys_dir / "ions.mdp", _minim_mdp(cfg), dry_run=dry_run)
        _run([gmx, "grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "hybrid.top",
              "-o", "ions.tpr", "-maxwarn", "2"], cwd=sys_dir, dry_run=dry_run)
        _run([gmx, "genion", "-s", "ions.tpr", "-o", "ions.gro", "-p", "hybrid.top",
              "-pname", "NA", "-nname", "CL", "-neutral",
              "-conc", str(scfg["ion_conc_M"]), "-seed", str(seed)],
             cwd=sys_dir, stdin="SOL\n", dry_run=dry_run)

        # 5. energy minimisation -- also relaxes the freed thiols, which start ~2 A apart
        #    because the crystal disulfide was broken rather than re-modelled.
        _write_mdp(sys_dir / "em.mdp", _minim_mdp(cfg), dry_run=dry_run)
        _run([gmx, "grompp", "-f", "em.mdp", "-c", "ions.gro", "-p", "hybrid.top",
              "-o", "em.tpr", "-maxwarn", "2"], cwd=sys_dir, dry_run=dry_run)
        _run([*mdrun_argv(cfg, gpu=False), "-deffnm", "em"], cwd=sys_dir, dry_run=dry_run,
             env={"OMP_NUM_THREADS": str(cfg["cluster"]["mdrun_ntomp"])})

        if not dry_run:
            (sys_dir / _DONE_MARKER).write_text(
                f"variant={variant} leg={leg} ff={ff} seed={seed}\n")
        _log(f"system ready: {sys_dir}")
        return sys_dir


# --------------------------------------------------------------------------- #
# .mdp generation                                                              #
# --------------------------------------------------------------------------- #
def _write_mdp(path: Path, body: str, dry_run: bool = False) -> Path:
    if not dry_run:
        path.write_text(body)
    return path


def _free_energy_lines(cfg: dict, window: int, nstdhdl: int | None = None) -> list[str]:
    """The ``free-energy`` + soft-core block, shared by minimisation and production.

    Minimisation needs this block too. Without ``free-energy = yes`` grompp builds the
    system from A-state parameters ALONE, so the atoms that appear only in state B are
    dummies -- no charge, no LJ -- and minimisation happily relaxes real atoms on top of
    them. Switching them on at a high lambda then explodes the system on step 1. See the
    "lambda-aware minimisation" note in CLAUDE.md.

    Args:
        cfg: pipeline config.
        window: lambda state this ``.mdp`` runs at (``init-lambda-state``).
        nstdhdl: dH/dl output stride. ``None`` for minimisation, which writes no dhdl.
    """
    fcfg = cfg["fep"]
    n_states = int(fcfg["lambda_windows"])
    lambdas = " ".join(f"{i / (n_states - 1):.4f}" for i in range(n_states))
    lines = [
        "free-energy              = yes",
        f"init-lambda-state        = {window}",
        f"fep-lambdas              = {lambdas}",
    ]
    if nstdhdl is not None:
        lines += [
            "calc-lambda-neighbors    = -1",   # energies at ALL states -> u_kn for MBAR
            f"nstdhdl                  = {nstdhdl}",
            "dhdl-print-energy        = total",
        ]
    # The two soft-core forms take DIFFERENT parameters. Beutler uses sc-alpha/power/sigma;
    # gapsys uses the sc-gapsys-* set and ignores (or rejects) the Beutler ones. Emitting
    # both is a grompp error, so pick one.
    sc = str(fcfg["sc_function"]).lower()
    if sc == "gapsys":
        lines += [
            "sc-function              = gapsys",
            f"sc-gapsys-scale-linpoint-lj = {fcfg['sc_gapsys_scale_linpoint_lj']}",
            f"sc-gapsys-scale-linpoint-q  = {fcfg['sc_gapsys_scale_linpoint_q']}",
            f"sc-gapsys-sigma-lj          = {fcfg['sc_gapsys_sigma_lj']}",
        ]
    elif sc == "beutler":
        lines += [
            "sc-function              = beutler",
            f"sc-alpha                 = {fcfg['sc_alpha']}",
            f"sc-power                 = {fcfg['sc_power']}",
            f"sc-sigma                 = {fcfg['sc_sigma']}",
            # sc-coul defaults to NO. One fep-lambdas vector drives charges and LJ
            # together, so without this an appearing atom carries a linearly-scaled
            # partial charge while its LJ core is already softened -- a neighbour can
            # penetrate the soft core and feel a real charge. Near-singular electrostatics
            # at intermediate lambda, worst in a packed core. pmx's published protein
            # protocol sets it; omitting it was a silent deviation.
            "sc-coul                  = yes",
        ]
    else:
        raise ValueError(f"fep.sc_function={sc!r} -- expected 'beutler' or 'gapsys'")
    return lines


def _minim_mdp(cfg: dict, window: int | None = None) -> str:
    """Steepest-descent minimisation.

    Args:
        cfg: pipeline config.
        window: when given, minimise AT THIS LAMBDA STATE rather than at the A state.
            The shared system build passes ``None`` (solvation/ions only care about the
            physical topology); every production window passes its own index.
    """
    fcfg = cfg["fep"]
    lines = [
        "integrator               = steep",
        "nsteps                   = 5000",
        "emtol                    = 100",
        f"coulombtype              = {fcfg['nonbonded_method']}",
        f"rvdw                     = {fcfg['nonbonded_cutoff_nm']}",
        f"rcoulomb                 = {fcfg['nonbonded_cutoff_nm']}",
    ]
    if window is not None:
        lines += ["", "; --- free energy: minimise at THIS window's lambda ---"]
        lines += _free_energy_lines(cfg, window)
    lines.append("")
    return "\n".join(lines)


def production_mdp(cfg: dict, window: int, seed: int, smoke: bool = False) -> str:
    """The production ``.mdp`` for one lambda window.

    A single ``fep-lambdas`` vector drives all interaction types, so ``init-lambda-state``
    alone selects this window. ``nstdhdl`` is set so the number of dH/dl records matches
    ``fep.frames_per_window``, and ``calc-lambda-neighbors = -1`` writes the energy at
    EVERY state each time -- that is what makes the output an MBAR ``u_kn`` rather than
    just a TI gradient.
    """
    fcfg = cfg["fep"]
    dt_ps = float(fcfg["timestep_fs"]) / 1000.0
    if smoke:
        production_steps, frames = 500, 10
        equil_steps = 100
    else:
        production_steps = int(round(float(fcfg["ns_per_window"]) * 1000.0 / dt_ps))
        frames = int(fcfg["frames_per_window"])
        equil_steps = int(round(float(fcfg["equilibration_ns"]) * 1000.0 / dt_ps))
    nsteps = production_steps + equil_steps
    nstdhdl = max(1, production_steps // frames)

    lines = [
        "; production FEP window -- generated by src/fep/pmx_engine.py",
        "integrator               = sd",          # stochastic dynamics = Langevin thermostat
        f"dt                       = {dt_ps}",
        f"nsteps                   = {nsteps}",
        f"nstcalcenergy            = {nstdhdl}",
        "nstlog                   = 5000",
        f"nstenergy                = {nstdhdl}",
        "nstxout-compressed       = 0",           # cluster.trajectory_retention: estimates only
        "",
        "; --- neighbour search / electrostatics ---",
        "cutoff-scheme            = Verlet",
        f"coulombtype              = {fcfg['nonbonded_method']}",
        f"rcoulomb                 = {fcfg['nonbonded_cutoff_nm']}",
        f"rvdw                     = {fcfg['nonbonded_cutoff_nm']}",
        "DispCorr                 = EnerPres",
        "",
        "; --- temperature / pressure ---",
        "tc-grps                  = System",
        f"tau-t                    = {1.0 / float(fcfg['collision_rate_per_ps']):.3f}",
        f"ref-t                    = {fcfg['temperature_K']}",
        f"gen-vel                  = yes",
        f"gen-temp                 = {fcfg['temperature_K']}",
        f"gen-seed                 = {seed}",
        f"ld-seed                  = {seed}",
    ]
    if fcfg["ensemble"] == "NPT":
        lines += [
            "pcoupl                   = C-rescale",
            "pcoupltype               = isotropic",
            "tau-p                    = 1.0",
            f"ref-p                    = {fcfg['pressure_atm']}",
            "compressibility          = 4.5e-5",
        ]
    lines += [
        "",
        "; --- constraints (>2 fs additionally needs HMR via pdb2gmx -heavyh; see fep.heavyh) ---",
        f"constraints              = {'h-bonds' if fcfg['constraints'] == 'HBonds' else 'none'}",
        "constraint-algorithm     = lincs",
        "lincs-order              = 6",
        "",
        "; --- free energy ---",
    ]
    lines += _free_energy_lines(cfg, window, nstdhdl)
    # NOTE: no couple-moltype here (we perturb a residue in place, we do not decouple a
    # whole molecule), so couple-intramol would be meaningless -- grompp warns on it.
    lines += [
        "",
        f"; equilibration steps discarded downstream: {equil_steps}",
        "",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# dhdl.xvg -> u_kn                                                             #
# --------------------------------------------------------------------------- #
def dhdl_to_u_kn(xvg: Path, temperature_K: float, n_states: int) -> np.ndarray:
    """Convert a GROMACS ``dhdl.xvg`` into ``(n_states, n_samples)`` reduced potentials.

    alchemlyb returns a DataFrame indexed by (time, lambda) whose columns are the lambda
    states, already in reduced (kT) units -- exactly MBAR's u_kn, transposed.

    Raises:
        ValueError: If the file yields a different number of states than configured,
            which means these windows are not the same ladder.
    """
    from alchemlyb.parsing.gmx import extract_u_nk

    u_nk = extract_u_nk(str(xvg), T=temperature_K)
    u_kn = np.asarray(u_nk.values, dtype=float).T     # (states, samples)
    if u_kn.shape[0] != n_states:
        raise ValueError(
            f"{xvg}: dhdl has {u_kn.shape[0]} lambda states, config says {n_states}. "
            "fep-lambdas in the .mdp and fep.lambda_windows disagree."
        )
    return u_kn


# --------------------------------------------------------------------------- #
# One window                                                                   #
# --------------------------------------------------------------------------- #
def run_pmx_window(cfg: dict, variant: str, leg: str, window: int, rep: int,
                   smoke: bool = False, dry_run: bool = False) -> dict:
    """Run one lambda window under GROMACS and return reduced potentials at every state.

    Returns the ``src.fep.window`` NPZ dict: ``lambda_index``, ``n_states``,
    ``u_kn_window`` of shape ``(n_states, n_frames)``, and ``provenance``.
    """
    fcfg = cfg["fep"]
    n_states = int(fcfg["lambda_windows"])
    gmx = gmx_command(cfg)

    _log(f"window start: {variant}/{leg} w{window} r{rep} smoke={smoke} n_states={n_states}")
    sys_dir = build_system(cfg, variant, leg, dry_run=dry_run)

    run_dir = ROOT / "results" / "fep" / variant / leg / f"w{window}_r{rep}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Sampling seed: distinct per replicate -> independent AND reproducible replicates.
    seed = stable_seed(variant, leg, window, rep, "sampling") % (2**31 - 1) or 1

    # Minimise AT THIS WINDOW'S LAMBDA before production. The system-wide em.gro was
    # minimised with free-energy off, i.e. at the A state, where the appearing atoms are
    # dummies -- real atoms get relaxed into them, and switching them on at a high lambda
    # detonates the system on step 1 (A4V folded w17 did exactly that, all replicates).
    # Cheap insurance: EM is seconds against ~18 min of production.
    em_gro = run_dir / "em.gro"
    if not em_gro.exists():
        _write_mdp(run_dir / "em.mdp", _minim_mdp(cfg, window), dry_run)
        _run([gmx, "grompp", "-f", "em.mdp",
              "-c", str(sys_dir / "em.gro"), "-p", str(sys_dir / "hybrid.top"),
              "-o", "em.tpr", "-maxwarn", "2"], cwd=run_dir, dry_run=dry_run)
        _run([*mdrun_argv(cfg, gpu=False), "-deffnm", "em"], cwd=run_dir, dry_run=dry_run,
             env={"OMP_NUM_THREADS": str(cfg["cluster"]["mdrun_ntomp"])})

    # Fingerprint the settings, not just the engine. `provenance` catches a mock or a
    # different engine; it does NOT catch the same engine run with sc-coul flipped or a
    # different nstdhdl. Mixing those inside one MBAR estimate is silent and wrong, which
    # is exactly the class of error the provenance rule exists to prevent. The seed is
    # excluded so replicates of one protocol still agree.
    mdp_text = production_mdp(cfg, window, seed, smoke)
    fingerprint = "\n".join(
        line for line in mdp_text.splitlines()
        if not line.startswith(("gen-seed", "ld-seed", "init-lambda-state", ";"))
    )
    protocol = hashlib.sha256(fingerprint.encode()).hexdigest()[:16]
    _write_mdp(run_dir / "prod.mdp", mdp_text, dry_run)

    _run([gmx, "grompp", "-f", "prod.mdp",
          "-c", str(em_gro if not dry_run else sys_dir / "em.gro"),
          "-p", str(sys_dir / "hybrid.top"),
          "-o", "prod.tpr", "-maxwarn", "2"], cwd=run_dir, dry_run=dry_run)

    # -cpi/-cpo: GROMACS resumes from its own checkpoint if the job was killed (CLAUDE.md).
    # -cpi/-append only make sense once a checkpoint exists; on a fresh run mdrun
    # errors out if they are passed without prod.cpt present.
    resume = ["-cpi", "prod.cpt", "-append"] if (run_dir / "prod.cpt").exists() else []
    # -cpt is in MINUTES and defaults to 15 -- about as long as a whole folded window, so
    # without it a killed window has no checkpoint to resume from (CLAUDE.md: assume
    # windows die). fep.checkpoint_interval_steps is the perses engine's key, not this one.
    cpt = ["-cpt", str(fcfg["checkpoint_interval_min"])]
    _run([*mdrun_argv(cfg), "-deffnm", "prod",
          *resume, "-cpo", "prod.cpt", *cpt,
          "-dhdl", "dhdl.xvg"], cwd=run_dir, dry_run=dry_run,
         env={"OMP_NUM_THREADS": str(cfg["cluster"]["mdrun_ntomp"])})

    if dry_run:
        return {"lambda_index": window, "n_states": n_states,
                "u_kn_window": np.zeros((n_states, 0)), "provenance": fcfg["framework"],
                "protocol": protocol}

    u_kn = dhdl_to_u_kn(run_dir / "dhdl.xvg", float(fcfg["temperature_K"]), n_states)
    dt_ps = float(fcfg["timestep_fs"]) / 1000.0
    if smoke:
        production_steps, frames = 500, 10
        equil_steps = 100
    else:
        production_steps = int(round(float(fcfg["ns_per_window"]) * 1000.0 / dt_ps))
        frames = int(fcfg["frames_per_window"])
        equil_steps = int(round(float(fcfg["equilibration_ns"]) * 1000.0 / dt_ps))
    nstdhdl = max(1, production_steps // frames)
    discard = (equil_steps + nstdhdl - 1) // nstdhdl
    u_kn = u_kn[:, discard:]
    if not np.all(np.isfinite(u_kn)):
        bad = int((~np.isfinite(u_kn)).sum())
        raise FloatingPointError(
            f"{variant}/{leg} w{window} r{rep}: {bad} non-finite reduced potentials "
            "(blown-up dynamics). Not writing a result."
        )
    _log(f"window done: {variant}/{leg} w{window} r{rep} u_kn{u_kn.shape} finite")
    return {"lambda_index": window, "n_states": n_states, "u_kn_window": u_kn,
            "provenance": fcfg["framework"], "protocol": protocol}


def main() -> None:
    """``python -m src.fep.pmx_engine --verify`` / ``--dry-run`` for first-light checks."""
    import argparse

    from src.prep.build import load_config

    p = argparse.ArgumentParser(description="GROMACS+pmx engine checks (Stage 3).")
    p.add_argument("--config", default=str(ROOT / "config" / "pipeline.yaml"))
    p.add_argument("--verify", action="store_true", help="check gmx/pmx/force field only")
    p.add_argument("--dry-run", action="store_true", help="print the command sequence")
    p.add_argument("--variant", default="A4V")
    p.add_argument("--leg", default="folded")
    p.add_argument("--window", type=int, default=0)
    p.add_argument("--rep", type=int, default=0)
    args = p.parse_args()

    cfg = load_config(args.config)
    if args.verify:
        for k, v in verify_tools(cfg).items():
            print(f"  {k:16} {v}")
        return
    run_pmx_window(cfg, args.variant, args.leg, args.window, args.rep,
                   smoke=True, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
