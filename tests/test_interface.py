"""Guards the geometric dimer-interface classifier (D4).

Two layers:
- A committed synthetic fixture (``tests/fixtures/mini_dimer.pdb``) pins the core
  logic deterministically, with no large-file or network dependency.
- A check against the real reference structure (3ECU) that SKIPS when the
  gitignored ``data/structures/3ECU.pdb`` is absent -- it asserts a canonical SOD1
  interface residue resolves to ``dimer`` and a buried-core residue to ``monomer``,
  so the rule stays reproducible (CLAUDE.md rule 3).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from src.panel.interface import (
    classify_oligomer,
    detect_partner_chain,
    interface_residues,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "mini_dimer.pdb"
CONFIG_PATH = ROOT / "config" / "pipeline.yaml"
STRUCTURE_3ECU = ROOT / "data" / "structures" / "3ECU.pdb"


@pytest.fixture(scope="module")
def structure_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)["structure"]


# ---- synthetic fixture: deterministic core logic ---------------------------

def test_fixture_partner_detection():
    assert detect_partner_chain(str(FIXTURE), "A", 5.0) == "B"


def test_fixture_interface_set():
    # residue 1 is placed within 5 A of chain B; residue 2 is 50 A away.
    assert interface_residues(str(FIXTURE), 5.0, ref_chain="A") == {1}


def test_fixture_classification():
    assert classify_oligomer(1, str(FIXTURE), 5.0, default="monomer") == "dimer"
    assert classify_oligomer(2, str(FIXTURE), 5.0, default="monomer") == "monomer"


def test_tight_cutoff_excludes_all():
    # nearest inter-chain contact is 2 A; a 1 A cutoff leaves no interface residues.
    assert interface_residues(str(FIXTURE), 1.0, ref_chain="A") == set()


# ---- real reference structure (skips if the file is not present) -----------

@pytest.mark.skipif(
    not STRUCTURE_3ECU.exists(),
    reason="data/structures/3ECU.pdb not fetched (gitignored)",
)
def test_3ecu_canonical_interface(structure_cfg):
    cutoff = structure_cfg["interface_cutoff_A"]
    default = structure_cfg["default_oligomer"]
    pdb = str(STRUCTURE_3ECU)

    iface = interface_residues(pdb, cutoff, ref_chain="A")
    # sanity: within the mature range and non-trivial
    assert iface and all(1 <= r <= 153 for r in iface)

    # canonical SOD1 dimer-interface residues (50s loop + C-terminal strand)
    assert classify_oligomer(51, pdb, cutoff, default=default) == "dimer"
    assert classify_oligomer(153, pdb, cutoff, default=default) == "dimer"

    # A4V and G93A sit in the buried beta-barrel core, not the interface
    assert classify_oligomer(4, pdb, cutoff, default=default) == "monomer"
    assert classify_oligomer(93, pdb, cutoff, default=default) == "monomer"
