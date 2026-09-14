"""CPU-only checks for the forensic reference sensitivity analysis."""

import numpy as np
import pytest

from src.analysis.reference_sensitivity import sensitivity_metrics


PRED = np.asarray([
    0.8292101987940607,
    2.5907683805493433,
    3.5441177106402493,
    1.2586774585329301,
    1.20360944032974,
    4.981916129859519,
    2.8309478830758548,
])


def test_gate_of_record_metrics_are_preserved():
    obs = np.asarray([0.37, 1.25, 1.62, 2.43, 3.70, 4.05, 7.00])
    metrics = sensitivity_metrics(PRED, obs)
    assert metrics["pearson"] == pytest.approx(0.32601995)
    assert metrics["rmse"] == pytest.approx(2.12348058)


def test_normalization_and_temperature_are_separable():
    normalized = np.asarray([0.37, 1.25, 1.62, 2.43, 1.85, 4.05, 3.50])
    temp_only = np.asarray([0.37, 1.25, 1.62, 2.43, 2.00, 4.05, 5.10])
    both = np.asarray([0.37, 1.25, 1.62, 2.43, 1.00, 4.05, 2.55])
    assert sensitivity_metrics(PRED, normalized)["pearson"] == pytest.approx(0.65696621)
    assert sensitivity_metrics(PRED, temp_only)["pearson"] == pytest.approx(0.53931562)
    assert sensitivity_metrics(PRED, both)["pearson"] == pytest.approx(0.75426141)
    assert sensitivity_metrics(PRED, both)["rmse"] == pytest.approx(1.07381932)
