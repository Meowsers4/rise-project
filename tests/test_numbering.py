"""Guards the mature<->precursor numbering offset (CLAUDE.md rule 3).

The offset must live ONLY in ``config/pipeline.yaml`` -> ``project.mature_offset``.
These tests load it from there and check round-trips + known variants, so a stray
hardcoded index or config drift is caught by ``pytest tests/``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from src.panel.numbering import (
    mature_to_precursor,
    parse_variant,
    precursor_to_mature,
)

CONFIG_PATH = Path(__file__).resolve().parents[1] / "config" / "pipeline.yaml"


@pytest.fixture(scope="module")
def project_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)["project"]


def test_offset_is_a_config_int(project_cfg):
    assert isinstance(project_cfg["mature_offset"], int)
    assert isinstance(project_cfg["n_residues_mature"], int)


def test_parse_known_variants():
    assert parse_variant("A4V") == ("A", 4, "V")
    assert parse_variant("G93A") == ("G", 93, "A")
    assert parse_variant(" G37R ") == ("G", 37, "R")  # whitespace tolerated


@pytest.mark.parametrize("bad", ["A4", "4V", "AV4", "B4V", "A4Z", "a4v", "", "A0V"])
def test_parse_rejects_malformed(bad):
    with pytest.raises(ValueError):
        parse_variant(bad)


def test_precursor_offset_applied(project_cfg):
    offset = project_cfg["mature_offset"]
    _, pos, _ = parse_variant("A4V")  # mature residue 4
    assert mature_to_precursor(pos, offset) == 4 + offset


def test_roundtrip_is_identity(project_cfg):
    offset = project_cfg["mature_offset"]
    for mature_pos in (1, 4, 93, project_cfg["n_residues_mature"]):
        precursor = mature_to_precursor(mature_pos, offset)
        assert precursor_to_mature(precursor, offset) == mature_pos


def test_control_positions_within_mature_length(project_cfg):
    n = project_cfg["n_residues_mature"]
    for name in ("A4V", "G93A", "G37R"):
        _, pos, _ = parse_variant(name)
        assert 1 <= pos <= n
