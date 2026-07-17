"""Stage 2 prescreen: parsers, mock-backed per-variant run, controls correlation.

No FoldX/Rosetta binaries required -- the licensed calls are behind a backend the
tests replace with MockBackend / direct construction.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.prescreen.run import (
    MockBackend,
    correlate_vs_experimental,
    parse_foldx_ddg,
    parse_rosetta_cartesian,
    prescreen_variant,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS_CSV = ROOT / "data" / "variants.csv"

CFG = {
    "panel": {"csv": "data/variants.csv"},
    "prescreen": {
        "methods": ["foldx", "rosetta_cartesian_ddg"],
        "runs_per_variant": 5,
        "rosetta_reu_per_kcal": 2.94,
    },
}

_FOLDX_DIF = (
    "FoldX 5.0 (c) ...\n"
    "\n"
    "Pdb\ttotal energy\tBackbone Hbond\t...\n"
    "A4V_1.pdb\t1.85\t0.0\t0.1\n"
)

_ROSETTA_DDG = (
    "COMPLEX:   Round1:  WT_A4V     -520.10   fa_atr -100 ...\n"
    "COMPLEX:   Round1:  MUT_A4V    -515.20   fa_atr -100 ...\n"
    "COMPLEX:   Round2:  WT_A4V     -519.90   fa_atr -100 ...\n"
    "COMPLEX:   Round2:  MUT_A4V    -514.80   fa_atr -100 ...\n"
)


def test_parse_foldx_ddg():
    assert parse_foldx_ddg(_FOLDX_DIF) == 1.85


def test_parse_foldx_ddg_errors_when_no_data():
    with pytest.raises(ValueError):
        parse_foldx_ddg("header only\nno data rows\n")


def test_parse_rosetta_cartesian_sign_and_scaling():
    # mean(MUT) - mean(WT) = -515.0 - (-520.0) = +5.0 REU; /2.94 -> +1.70 kcal, destabilizing
    ddg = parse_rosetta_cartesian(_ROSETTA_DDG, reu_per_kcal=2.94)
    assert ddg == pytest.approx(5.0 / 2.94, abs=1e-6)
    assert ddg > 0  # positive = destabilizing


def test_parse_rosetta_missing_rows_errors():
    with pytest.raises(ValueError):
        parse_rosetta_cartesian("COMPLEX: Round1: WT_x -1.0\n", reu_per_kcal=2.94)


@pytest.mark.skipif(not VARIANTS_CSV.exists(), reason="variants.csv not built")
def test_prescreen_variant_writes_json(tmp_path):
    out = tmp_path / "A4V.json"
    result = prescreen_variant(CFG, "A4V", out, backend=MockBackend())
    assert out.exists()
    saved = json.loads(out.read_text())
    assert saved["variant"] == "A4V"
    assert set(saved["methods"]) == {"foldx", "rosetta_cartesian_ddg"}
    assert saved["methods"]["foldx"]["n"] == 5      # runs_per_variant averaging
    assert "consensus_ddg" in saved


def test_correlate_vs_experimental(tmp_path):
    # Build a tiny controls csv and perfectly-correlated predictions.
    csv_path = tmp_path / "variants.csv"
    csv_path.write_text(
        "variant,bucket,exp_ddg\n"
        "A4V,positive_control,1.62\n"
        "G93A,positive_control,2.43\n"
        "G93V,positive_control,7.00\n"
    )
    results = [
        {"variant": "A4V", "consensus_ddg": 1.6, "methods": {}},
        {"variant": "G93A", "consensus_ddg": 2.5, "methods": {}},
        {"variant": "G93V", "consensus_ddg": 7.1, "methods": {}},
    ]
    stats = correlate_vs_experimental(results, csv_path)
    assert stats["n"] == 3
    assert stats["pearson"] > 0.99  # monotone, near-perfect


def test_correlate_too_few_points():
    stats = correlate_vs_experimental([{"variant": "A4V", "consensus_ddg": 1.0}], VARIANTS_CSV)
    # A4V is a real control, but n<3 overall -> None
    assert stats["pearson"] is None
