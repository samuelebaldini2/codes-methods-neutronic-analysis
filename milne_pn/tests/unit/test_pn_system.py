"""
Unit tests for milne_pn.algorithms.pn_system.PNSystem.
"""

import numpy as np
import pytest

from milne_pn.algorithms.pn_system import PNSystem


def test_rejects_even_order():
    with pytest.raises(ValueError):
        PNSystem(4)


def test_rejects_nonpositive_order():
    with pytest.raises(ValueError):
        PNSystem(-1)


def test_matrix_shapes():
    N = 7
    pn = PNSystem(N)
    assert pn.A.shape == (N + 1, N + 1)
    assert pn.D.shape == (N + 1, N + 1)
    assert pn.H.shape == (N + 1, N + 1)


def test_A_is_zero_diagonal_tridiagonal():
    pn = PNSystem(5)
    A = pn.A
    n = A.shape[0]
    assert np.allclose(np.diag(A), 0.0)
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1:
                assert A[i, j] == 0.0


def test_A_known_entries():
    pn = PNSystem(3)
    A = pn.A
    assert A[0, 1] == pytest.approx(1.0)
    assert A[1, 0] == pytest.approx(1 / 3)
    assert A[1, 2] == pytest.approx(2 / 3)
    assert A[2, 1] == pytest.approx(2 / 5)
    assert A[2, 3] == pytest.approx(3 / 5)
    assert A[3, 2] == pytest.approx(3 / 7)
    assert A[3, 3] == 0.0


def test_D_diagonal():
    pn = PNSystem(5)
    D = pn.D
    assert D[0, 0] == 0.0
    for l in range(1, pn.size):
        assert D[l, l] == 1.0
    assert np.count_nonzero(D - np.diag(np.diag(D))) == 0


def test_H_matches_known_low_order_integrals():
    pn = PNSystem(3)
    H = pn.H
    assert H[0, 0] == pytest.approx(1.0)
    assert H[1, 0] == pytest.approx(0.5)
    assert H[1, 1] == pytest.approx(1 / 3)
    assert np.allclose(H, H.T, atol=1e-12)


def test_H_positive_definite():
    pn = PNSystem(9)
    eigvals = np.linalg.eigvalsh(pn.H)
    assert np.all(eigvals > 0)
