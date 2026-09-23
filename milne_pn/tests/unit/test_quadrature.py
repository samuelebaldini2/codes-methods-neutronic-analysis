"""
Unit tests for milne_pn.utils.quadrature.
"""

import numpy as np
import pytest

from milne_pn.utils.quadrature import gauss_legendre_half_range, legendre_matrix


def test_half_range_weights_sum_to_one():
    mu, wq = gauss_legendre_half_range(20)
    assert np.sum(wq) == pytest.approx(1.0)


def test_half_range_nodes_within_unit_interval():
    mu, wq = gauss_legendre_half_range(20)
    assert np.all(mu >= 0.0) and np.all(mu <= 1.0)


def test_half_range_integrates_polynomials_exactly():
    n = 10
    mu, wq = gauss_legendre_half_range(n)
    for power in range(2 * n - 1):
        approx = np.sum(wq * mu ** power)
        exact = 1.0 / (power + 1)
        assert approx == pytest.approx(exact, rel=1e-10)


def test_legendre_matrix_endpoint_values():
    mu = np.array([1.0, -1.0, 0.0])
    P = legendre_matrix(4, mu)
    assert np.allclose(P[:, 0], 1.0)
    assert np.allclose(P[:, 1], [(-1) ** l for l in range(5)])
    assert P[0, 2] == pytest.approx(1.0)
    assert P[1, 2] == pytest.approx(0.0)
    assert P[2, 2] == pytest.approx(-0.5)


def test_legendre_matrix_orthogonality_full_range():
    n = 30
    x, w = np.polynomial.legendre.leggauss(n)
    P = legendre_matrix(6, x)
    for k in range(7):
        for l in range(7):
            val = np.sum(P[k] * P[l] * w)
            expected = 2.0 / (2 * k + 1) if k == l else 0.0
            assert val == pytest.approx(expected, abs=1e-10)
