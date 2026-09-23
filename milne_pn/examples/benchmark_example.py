#!/usr/bin/env python3
"""
examples/benchmark_example.py
--------------------------------
Runs the convergence benchmark across PN order, and the Runge/Gibbs
phenomenon analysis, then renders all four plots into examples/output/.
"""

import os

import numpy as np

from milne_pn.benchmarks.benchmark import convergence_benchmark
from milne_pn.algorithms.milne_solver import MilneSolver
from milne_pn.visualization import (
    plot_convergence,
    plot_flux_profile,
    plot_gibbs_ringing,
    plot_runge_phenomenon,
)

Z0_EXACT = 0.710446


def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)

    result = convergence_benchmark(n_max=59, z0_exact=Z0_EXACT)
    plot_convergence(result["N"], result["z0"], Z0_EXACT,
                      os.path.join(out_dir, "convergence.png"))

    N_plot = 31
    solver = MilneSolver(N_plot)
    x = np.linspace(0, 6, 300)
    plot_flux_profile(x, solver.phi0(x), solver.z0,
                       os.path.join(out_dir, "flux_profile.png"), N_label=N_plot)

    plot_gibbs_ringing([3, 9, 21, 41], os.path.join(out_dir, "gibbs_ringing.png"))
    plot_runge_phenomenon([5, 10, 15, 20], os.path.join(out_dir, "runge_phenomenon.png"))

    n_last = result["N"][-1]
    z0_last = result["z0"][-1]
    print(f"Wrote 4 plots to {out_dir}/")
    print(f"z0(N={n_last}) = {z0_last:.6f}  (error {z0_last - Z0_EXACT:+.2e})")


if __name__ == "__main__":
    main()
