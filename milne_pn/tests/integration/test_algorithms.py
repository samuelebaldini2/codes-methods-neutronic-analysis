"""
Integration test: the full pipeline (PNSystem -> MilneSolver -> z0), and the
Runge/Gibbs-phenomenon analysis built on top of it, exercised end to end.
"""

import numpy as np
import pytest

from milne_pn.algorithms.milne_solver import MilneSolver
from milne_pn.algorithms.runge_analysis import (
    reconstruct_boundary_flux,
    gibbs_overshoot,
    runge_interpolation_demo,
)
from milne_pn.benchmarks.benchmark import convergence_benchmark

Z0_EXACT = 0.710446


def test_convergence_benchmark_end_to_end():
    result = convergence_benchmark(n_max=41, z0_exact=Z0_EXACT)
    assert result["N"][0] == 1
    assert result["N"][-1] == 41
    assert np.all(np.diff(result["z0"]) > 0)
    assert abs(result["z0"][-1] - Z0_EXACT) < 1e-3


def test_high_order_matches_case_to_four_digits():
    solver = MilneSolver(79)
    assert solver.z0 == pytest.approx(Z0_EXACT, abs=2e-5)


def test_gibbs_overshoot_pipeline_runs_and_is_bounded():
    for N in (5, 15, 25):
        solver = MilneSolver(N)
        mu, psi0 = reconstruct_boundary_flux(solver, n_mu=500)
        assert mu.shape == psi0.shape
        stats = gibbs_overshoot(solver, n_mu=500)
        assert stats["N"] == N
        assert 0.0 <= stats["relative_overshoot"] < 1.0


def test_gibbs_overshoot_does_not_vanish_with_order():
    orders = [9, 21, 41, 61]
    rel = [gibbs_overshoot(MilneSolver(N), n_mu=1000)["relative_overshoot"] for N in orders]
    assert rel[-1] > 0.5 * rel[0]


def test_runge_interpolation_demo_shows_growing_error():
    result = runge_interpolation_demo(degrees=(5, 10, 15, 20), n_eval=500)
    errs = [result["max_abs_error"][d] for d in (5, 10, 15, 20)]
    assert errs[-1] > errs[0]
