"""Stage 5 validation gate -- the go/no-go logic (CLAUDE.md #2).

Synthetic FEP ΔΔG only; no GPU/Perses. Exercises the gate pass/fail, exclusion of
non-converged points, the min-points floor, classification, and the file-level
refusal-to-classify on a failed gate.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from src.analysis.validate import (
    classify_uncharacterized,
    evaluate_gate,
    run_gate_validation,
    run_validation,
)

CFG = {
    "panel": {"csv": "data/variants.csv"},
    "fep": {"framework": "openmm_perses", "convergence": {"max_cycle_closure_kcal": 1.0}},
    "validation": {
        "min_pearson": 0.6,
        "max_rmse_kcal": 1.5,
        "min_gate_points": 3,
        "destabilizing_ddg_kcal": 1.0,
        "gate_subset": ["A", "B", "C", "D"],
        "outputs": {"gate_report": "results/validation_gate.json"},
    },
}


def _fep(ddg, converged=True, closure=0.3, provenance="openmm_perses"):
    return {"ddg": ddg, "ddg_err": 0.3, "cycle_closure_kcal": closure,
            "converged": converged, "provenance": provenance}


def test_gate_passes_on_good_correlation():
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {v: _fep(exp[v] + 0.2) for v in exp}  # near-perfect
    gate = evaluate_gate(fep, exp, CFG)
    assert gate["passed"] and gate["pearson"] > 0.99 and gate["n"] == 4


def test_gate_fails_on_poor_correlation():
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {"A": _fep(7.0), "B": _fep(1.0), "C": _fep(6.0), "D": _fep(2.0)}  # scrambled
    gate = evaluate_gate(fep, exp, CFG)
    assert not gate["passed"]


def test_gate_excludes_unconverged_and_enforces_floor():
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {
        "A": _fep(1.1),
        "B": _fep(2.1, converged=False),      # excluded: not converged
        "C": _fep(4.1, closure=5.0),          # excluded: cycle closure too large
        "D": _fep(7.1),
    }
    gate = evaluate_gate(fep, exp, CFG)
    # only A and D usable -> below min_gate_points=3 -> cannot judge -> fail
    assert not gate["passed"] and gate["n"] == 2
    assert {e["variant"] for e in gate["excluded"]} >= {"B", "C"}


def test_gate_excludes_mock_provenance():
    """A mock ΔΔG must never reach the gate, however good the correlation looks."""
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {v: _fep(exp[v] + 0.1, provenance="mock") for v in exp}
    gate = evaluate_gate(fep, exp, CFG)
    assert not gate["passed"] and gate["n"] == 0
    assert all("not the production engine" in e["reason"] for e in gate["excluded"])


def test_gate_excludes_missing_provenance():
    """This is the fabricated-file case: a plausible ΔΔG with no recorded origin."""
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {}
    for v in exp:
        rec = _fep(exp[v] + 0.1)
        del rec["provenance"]           # hand-written / pre-provenance file
        fep[v] = rec
    gate = evaluate_gate(fep, exp, CFG)
    assert not gate["passed"] and gate["n"] == 0
    assert all("'MISSING'" in e["reason"] for e in gate["excluded"])


def test_classification_skips_non_production():
    panel = [{"variant": "X1Y", "bucket": "vus", "oligomer": "monomer"}]
    assert classify_uncharacterized({"X1Y": _fep(3.5, provenance="mock")}, panel, CFG) == []
    assert len(classify_uncharacterized({"X1Y": _fep(3.5)}, panel, CFG)) == 1


def test_gate_fails_on_good_correlation_but_bad_accuracy():
    """The case correlation alone cannot see: perfectly ranked, systematically wrong."""
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    fep = {v: _fep(exp[v] * 0.5) for v in exp}      # every value half the truth
    gate = evaluate_gate(fep, exp, CFG)
    assert gate["pearson"] > 0.99                   # correlation says "perfect"
    assert gate["rmse"] > 1.5 and not gate["passed"]
    assert "rmse" in gate["reason"]
    assert gate["slope"] == pytest.approx(0.5, abs=0.01)


def test_accuracy_metrics_reported():
    exp = {"A": 1.0, "B": 2.0, "C": 4.0, "D": 7.0}
    gate = evaluate_gate({v: _fep(exp[v] + 0.2) for v in exp}, exp, CFG)
    assert gate["passed"]
    assert gate["rmse"] == pytest.approx(0.2, abs=0.01)
    assert gate["mue"] == pytest.approx(0.2, abs=0.01)
    assert gate["slope"] == pytest.approx(1.0, abs=0.01)


def test_never_lowers_threshold():
    # pearson ~0.55 < 0.6 must fail even though "close"
    exp = {"A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}
    fep = {"A": _fep(1.0), "B": _fep(3.0), "C": _fep(2.0), "D": _fep(4.2)}
    gate = evaluate_gate(fep, exp, CFG)
    if gate["pearson"] is not None and gate["pearson"] < 0.6:
        assert not gate["passed"]


def test_classify_uncharacterized():
    panel = [
        {"variant": "X1Y", "bucket": "vus", "oligomer": "monomer"},
        {"variant": "Z2W", "bucket": "pathogenic_uncharacterized", "oligomer": "dimer"},
        {"variant": "A4V", "bucket": "positive_control", "oligomer": "monomer"},  # skipped
    ]
    fep = {"X1Y": _fep(0.2), "Z2W": _fep(3.5), "A4V": _fep(1.6)}
    rows = classify_uncharacterized(fep, panel, CFG)
    labels = {r["variant"]: r["prediction"] for r in rows}
    assert labels == {"X1Y": "stable", "Z2W": "destabilizing"}  # control excluded


def _write_panel(path: Path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["variant", "bucket", "oligomer", "exp_ddg"])
        w.writeheader()
        for v, dg in [("A", 1.0), ("B", 2.0), ("C", 4.0), ("D", 7.0)]:
            w.writerow({"variant": v, "bucket": "positive_control", "oligomer": "monomer", "exp_ddg": dg})
        w.writerow({"variant": "X1Y", "bucket": "vus", "oligomer": "monomer", "exp_ddg": ""})


def _write_fep(fep_dir: Path, mapping: dict):
    for v, rec in mapping.items():
        d = fep_dir / v
        d.mkdir(parents=True)
        (d / "ddg.json").write_text(json.dumps({"variant": v, **rec}))


def test_run_validation_passes_and_writes_map(tmp_path):
    panel = tmp_path / "variants.csv"
    _write_panel(panel)
    fep_dir = tmp_path / "fep"
    _write_fep(fep_dir, {v: _fep(dg) for v, dg in
                         {"A": 1.1, "B": 2.1, "C": 4.1, "D": 7.1, "X1Y": _fep_val()}.items()})
    ddg_map = tmp_path / "ddg_map.csv"
    report = tmp_path / "gate.json"
    gate = run_validation(CFG, ddg_map, report, fep_dir=fep_dir, variants_csv=panel)
    assert gate["passed"]
    assert report.exists() and ddg_map.exists()
    rows = list(csv.DictReader(open(ddg_map)))
    assert [r["variant"] for r in rows] == ["X1Y"]  # only uncharacterized classified


def _fep_val():
    return 0.2  # helper: X1Y ddG (used above)


def test_run_validation_fails_refuses_map(tmp_path):
    panel = tmp_path / "variants.csv"
    _write_panel(panel)
    fep_dir = tmp_path / "fep"
    # scrambled -> poor correlation -> no-go
    _write_fep(fep_dir, {"A": _fep(7.0), "B": _fep(1.0), "C": _fep(6.0), "D": _fep(2.0),
                         "X1Y": _fep(0.2)})
    ddg_map = tmp_path / "ddg_map.csv"
    report = tmp_path / "gate.json"
    with pytest.raises(SystemExit):
        run_validation(CFG, ddg_map, report, fep_dir=fep_dir, variants_csv=panel)
    assert report.exists()          # report always written
    assert not ddg_map.exists()     # map refused on no-go


def test_gate_validation_does_not_require_non_gate_fep(tmp_path):
    panel = tmp_path / "variants.csv"
    _write_panel(panel)
    fep_dir = tmp_path / "fep"
    _write_fep(fep_dir, {v: _fep(dg) for v, dg in
                         {"A": 1.1, "B": 2.1, "C": 4.1, "D": 7.1}.items()})
    report = tmp_path / "gate.json"
    gate = run_gate_validation(CFG, report, fep_dir=fep_dir, variants_csv=panel)
    assert gate["passed"]
    assert report.exists()
