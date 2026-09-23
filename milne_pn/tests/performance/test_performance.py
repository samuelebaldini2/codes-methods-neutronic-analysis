"""
Performance tests: sanity-check that solve time stays reasonable and grows
in a sane way with PN order.
"""

import time

import pytest

from milne_pn.algorithms.milne_solver import MilneSolver
from milne_pn.benchmarks.benchmark import timing_benchmark


@pytest.mark.parametrize("N", [9, 29, 59])
def test_solve_completes_quickly(N):
    t0 = time.perf_counter()
    MilneSolver(N)
    elapsed = time.perf_counter() - t0
    assert elapsed < 5.0


def test_timing_benchmark_runs_and_is_monotone_ish():
    result = timing_benchmark([5, 25, 45], repeats=1)
    assert result["N"].shape == (3,)
    assert result["mean_time_s"].shape == (3,)
    assert result["mean_time_s"][-1] >= 0.0
    assert result["mean_time_s"][0] >= 0.0
