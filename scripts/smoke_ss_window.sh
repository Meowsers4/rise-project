#!/bin/bash -l
#
# smoke_ss_window.sh -- one-off verification job for the G93A disulfide diagnostic.
# DIAGNOSTIC BRANCH ONLY (diag/g93a-ss); not part of the pipeline, not a Snakemake rule.
#
# Builds the G93A folded system under the current config and reports whether pdb2gmx
# formed the intended C57-C146 bond and ONLY that bond. On this branch
# fep.keep_disulfide_reduced is false, which disables assert_topology_disulfide_free --
# so this grep is the only verification that the topology is what we intend.
#
# Usage, from the repo root on the SCC:
#   mkdir -p logs/fep && qsub scripts/smoke_ss_window.sh
#   tail -f logs/fep/ss_smoke.o<jobid>
#
# It cannot be an array task: submit_array.sh's runtime guard requires the submitted
# array size to equal legs*windows*replicates (120), so `-t 1-1` exits 2 by design.
#
#$ -P rise-batteries
#$ -l gpus=1
#$ -l gpu_c=7.0
#$ -l gpu_type=L40S
#$ -l h_rt=2:00:00
#$ -pe omp 8
#$ -N ss_smoke
#$ -j y
#$ -cwd
#$ -o logs/fep/

# NOT -u: GMXRC reads unbound $shell/$GMXLDLIB (CLAUDE.md first-light trap 6).
set -eo pipefail

source scripts/scc_env.sh
echo "host=$(hostname) branch=$(git branch --show-current) CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
python -c "import yaml;c=yaml.safe_load(open('config/pipeline.yaml'))['fep'];print('keep_disulfide_reduced =',c['keep_disulfide_reduced'],'| ns_per_window =',c['ns_per_window'])"

OUT="${TMPDIR:-/tmp}/ss_smoke_w0_r0.npz"
python -m src.fep.window --variant G93A --leg folded --window 0 --rep 0 --smoke \
    --config config/pipeline.yaml --out "$OUT"

T=results/fep/G93A/folded/system_r0/hybrid.top
echo; echo "==================== TOPOLOGY VERDICT ===================="
if [[ ! -f "$T" ]]; then echo "FAIL: $T does not exist"; exit 1; fi

# pmx/GROMACS name a bridged cysteine CYS2 (amber) or CYX; a free thiol keeps HG.
BRIDGED=$(grep -cE '\bCYS2\b|\bCYX\b' "$T" || true)
echo "bridged-cysteine residue lines in hybrid.top: ${BRIDGED}"
echo "--- residue names for the four cysteines (6, 57, 111, 146) ---"
python - "$T" <<'PY'
import re, sys
# [ atoms ] lines: nr type resnr residue atom cgnr charge mass
names={}
for line in open(sys.argv[1]):
    p=line.split()
    if len(p)>=5 and p[0].isdigit() and p[2].isdigit():
        rn=int(p[2]); nm=p[3]
        if nm.upper().startswith(("CYS","CYX")): names.setdefault(rn,set()).add(nm)
cys=sorted(names)
print(f"cysteine-like residues found: {len(cys)}")
for rn in cys: print(f"  resnr {rn:4d}  {sorted(names[rn])}")
bridged=[rn for rn,ns in names.items() if any(n.upper() in ("CYS2","CYX") for n in ns)]
print(f"\nbridged: {sorted(bridged)}   free thiol: {sorted(set(cys)-set(bridged))}")
print("EXPECT exactly 2 bridged (the C57/C146 pair) and 2 free (C6, C111).")
print("VERDICT:", "PASS-shape" if len(bridged)==2 else f"CHECK -- {len(bridged)} bridged, not 2")
PY
echo "--- SS-bond directives, if the force field emits them ---"
grep -n -A5 "disulf\|SSBOND\|; *bonds" "$T" | head -20 || true
echo; echo "--- window written ---"
python -c "
import numpy as np,os; z=np.load(os.environ['OUT'])
print('shape', z['u_kn_window'].shape, '| protocol', str(z['protocol']), '| provenance', str(z['provenance']))
print('NOTE: the protocol hash is IDENTICAL to the 2SH baseline by design -- it cannot witness the redox state.')"
echo "=========================================================="

# --smoke stamps this run dir with a different protocol hash; the real array task for
# w0_r0 would then hit assert_resumable and refuse. Remove it, keep the system dir.
rm -rf results/fep/G93A/folded/w0_r0
echo "removed smoke run dir; system_r0 kept (identical under either protocol)"
