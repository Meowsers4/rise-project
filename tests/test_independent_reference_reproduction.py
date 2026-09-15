"""Tests for the clean-room reference-sensitivity reproduction."""

import pytest

from scripts.reproduce_reference_sensitivity_independent import metrics, ranks


PREDICTED = [
    0.8292101987940607,
    2.5907683805493433,
    3.5441177106402493,
    1.2586774585329301,
    1.20360944032974,
    4.981916129859519,
    2.8309478830758548,
]


def test_rank_calculation_handles_ties_independently():
    assert ranks([30.0, 10.0, 20.0, 20.0]) == [4.0, 1.0, 2.5, 2.5]


def test_gate_metrics_reproduce_documented_full_precision_values():
    observed = [0.37, 1.25, 1.62, 2.43, 3.70, 4.05, 7.00]
    result = metrics(PREDICTED, observed)
    assert result["pearson"] == pytest.approx(0.32601995)
    assert result["spearman"] == pytest.approx(0.5)
    assert result["rmse"] == pytest.approx(2.12348058)
    assert result["mue"] == pytest.approx(1.7846825197006642)


def test_combined_constant_dcp_metrics_reproduce_documented_values():
    observed = [0.37, 1.25, 1.62, 2.43, 1.00, 4.05, 2.55]
    result = metrics(PREDICTED, observed)
    assert result["pearson"] == pytest.approx(0.75426141)
    assert result["spearman"] == pytest.approx(0.82142857)
    assert result["rmse"] == pytest.approx(1.07381932)
    assert result["mue"] == pytest.approx(0.9016988978165481)
