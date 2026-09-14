"""Unit tests for the archived-input methods figure builders."""

from pathlib import Path

import pytest

from src.analysis.methods_figures import (
    disagreement_correlations,
    folded_disagreement_rows,
    g93a_redox_summary,
    inspect_cysteine_topology,
    parse_cys_topology,
)
from src.analysis.archive_manifest import aggregate_digest, sha256_file


def _record(variant, folded):
    return {
        "variant": variant,
        "legs": [
            {
                "leg": "folded",
                "replicate": index,
                "dg_kcal": dg,
                "hysteresis_kcal": hysteresis,
            }
            for index, (dg, hysteresis) in enumerate(folded)
        ],
    }


def test_disagreement_uses_mean_of_the_two_siblings():
    rows = folded_disagreement_rows({
        "X1A": _record("X1A", [(1.0, 0.1), (2.0, 0.2), (4.0, 0.3)]),
    })
    assert [row["disagreement_kcal"] for row in rows] == pytest.approx([2.0, 0.5, 2.5])


def test_correlations_are_computed_on_raw_diagnostics():
    rows = [
        {"hysteresis_kcal": 1.0, "disagreement_kcal": 4.0},
        {"hysteresis_kcal": 2.0, "disagreement_kcal": 3.0},
        {"hysteresis_kcal": 3.0, "disagreement_kcal": 2.0},
        {"hysteresis_kcal": 4.0, "disagreement_kcal": 1.0},
    ]
    pearson, spearman = disagreement_correlations(rows)
    assert pearson == pytest.approx(-1.0)
    assert spearman == pytest.approx(-1.0)


def _write_topology(path: Path, oxidized: bool) -> None:
    hg57 = "" if oxidized else "3 HS 57 CYS HG 3 0.1 1.0\n"
    hg146 = "" if oxidized else "5 HS 146 CYS HG 5 0.1 1.0\n"
    bond = "1 4 1\n" if oxidized else ""
    path.write_text(
        "[ atoms ]\n"
        "1 S 57 CYS SG 1 -0.1 32.0\n"
        f"{hg57}"
        "4 S 146 CYS SG 4 -0.1 32.0\n"
        f"{hg146}"
        "[ bonds ]\n"
        f"{bond}"
    )


def test_topology_parser_distinguishes_ss_and_2sh(tmp_path):
    reduced = tmp_path / "reduced.top"
    oxidized = tmp_path / "oxidized.top"
    _write_topology(reduced, oxidized=False)
    _write_topology(oxidized, oxidized=True)
    reduced_state = parse_cys_topology(reduced)
    oxidized_state = parse_cys_topology(oxidized)
    assert reduced_state["bonded"] is False
    assert all(reduced_state["hg_present"].values())
    assert oxidized_state["bonded"] is True
    assert not any(oxidized_state["hg_present"].values())
    assert inspect_cysteine_topology(oxidized)["sg_bonds"] == [(57, 146)]


def test_redox_summary_uses_folded_leg_only():
    baseline = _record("G93A", [(9.9, 0.1), (9.6, 0.1), (9.7, 0.1)])
    ss = _record("G93A", [(9.8, 0.1), (9.8, 0.1), (9.7, 0.1)])
    summary = g93a_redox_summary(baseline, ss)
    assert summary["shift"] == pytest.approx(0.033333333333333)


def test_archive_digest_is_order_independent_and_content_sensitive(tmp_path):
    source = tmp_path / "input.dat"
    source.write_bytes(b"archived input\n")
    checksum = sha256_file(source)
    first = {"root": "archive", "path": "a", "size_bytes": 15, "sha256": checksum}
    second = {"root": "archive", "path": "b", "size_bytes": 15, "sha256": checksum}
    assert aggregate_digest([first, second]) == aggregate_digest([second, first])
    changed = dict(second, sha256="0" * 64)
    assert aggregate_digest([first, second]) != aggregate_digest([first, changed])
