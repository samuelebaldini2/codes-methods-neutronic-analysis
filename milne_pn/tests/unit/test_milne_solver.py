"""
Unit tests for milne_pn.algorithms.milne_solver.MilneSolver.
"""

import numpy as np
import pytest

from milne_pn.algorithms.milne_solver import MilneSolver


def test_P1_matches_classical_diffusion_theory():
    solver = MilneSolver(1)
    assert solver.z0 == pytest.approx(2.0 / 3.0, abs=1e-8)
    assert solver.num_discrete_modes == 0


def test_spectrum_is_real():
    for N in (3, 9, 15, 25):
        solver = MilneSolver(N)
        assert solver.max_imag_residual < 1e-8


def test_number_of_discrete_modes():
    for N in (1, 3, 5, 9, 15, 21):
        solver = MilneSolver(N)
        assert solver.num_discrete_modes == (N - 1) // 2


def test_z0_converges_monotonically_upward_toward_case_value():
    z0_exact = 0.710446
    orders = [1, 3, 9, 21, 41]
    z0s = [MilneSolver(N).z0 for N in orders]
    for a, b in zip(z0s, z0s[1:]):
        assert b > a
        assert b <= z0_exact + 1e-6


def test_z0_within_known_tolerance_at_moderate_order():
    solver = MilneSolver(41)
    assert abs(solver.z0 - 0.710446) < 1e-3


def test_phi0_is_linear_far_from_boundary():
    solver = MilneSolver(21)
    x = np.array([20.0, 30.0, 40.0])
    phi0 = solver.phi0(x)
    expected = x + solver.z0
    assert np.allclose(phi0, expected, atol=1e-6)


def test_phi0_boundary_value_is_positive_and_less_than_asymptote():
    solver = MilneSolver(21)
    boundary_val = solver.phi0(0.0)[0]
    assert 0.0 < boundary_val < solver.z0


def test_marshak_conditions_are_satisfied():
    solver = MilneSolver(9)
    from milne_pn.utils.quadrature import gauss_legendre_half_range, legendre_matrix

    mu, wq = gauss_legendre_half_range(4 * solver.pn.size)
    psi0 = solver.angular_flux(0.0, mu)[0]
    P = legendre_matrix(solver.N, mu)

    for k in range(1, solver.N + 1, 2):
        moment = np.sum(P[k] * psi0 * wq)
        assert abs(moment) < 1e-8
