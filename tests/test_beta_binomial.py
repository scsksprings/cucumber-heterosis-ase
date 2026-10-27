#!/usr/bin/env python3
import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "pipeline/scripts/fit_ase_beta_binomial.py"
spec = importlib.util.spec_from_file_location("ase_model", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_single_group_direction_and_significance():
    mu, phi, pvalue, eta_se, converged = module.fit_single_group(
        np.array([90.0, 88.0, 92.0]), np.array([100.0, 100.0, 100.0])
    )
    assert converged
    assert mu > 0.85
    assert pvalue < 0.05
    assert np.isfinite(eta_se)


def test_single_group_balanced():
    mu, _, pvalue, _, converged = module.fit_single_group(
        np.array([49.0, 50.0, 51.0]), np.array([100.0, 100.0, 100.0])
    )
    assert converged
    assert abs(mu - 0.5) < 0.02
    assert pvalue > 0.05


def test_two_group_difference():
    y = np.array([75.0, 78.0, 76.0, 25.0, 22.0, 24.0])
    n = np.repeat(100.0, 6)
    group = np.array([0, 0, 0, 1, 1, 1])
    mu0, mu1, _, pvalue, converged = module.fit_two_groups(y, n, group)
    assert converged
    assert mu0 > 0.7 and mu1 < 0.3
    assert pvalue < 0.05
