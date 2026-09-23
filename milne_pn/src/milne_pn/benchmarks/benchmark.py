"""
benchmark.py
------------
Convergence and timing studies for the PN Milne solver, across a range of
odd orders N.
"""

import argparse
import time

import numpy as np

from ..algorithms.milne_solver import MilneSolver

Z0_EXACT = 0.710446


def convergence_benchmark(n_max=79, n_min=1, z0_exact=Z0_EXACT):
    if n_min % 2 == 0:
        n_min += 1
    if n_max % 2 == 0:
        n_max += 1

    Ns, z0s, errs, modes, max_im = [], [], [], [], []
    for N in range(n_min, n_max + 1, 2):
        solver = MilneSolver(N)
        Ns.append(N)
        z0s.append(solver.z0)
        errs.append(solver.z0 - z0_exact)
        modes.append(solver.num_discrete_modes)
        max_im.append(solver.max_imag_residual)

    return {
        "N": np.array(Ns),
        "z0": np.array(z0s),
        "error": np.array(errs),
        "num_discrete_modes": np.array(modes),
        "max_imag_residual": np.array(max_im),
        "z0_exact": z0_exact,
    }


def timing_benchmark(orders, repeats=3):
    Ns, mean_times, std_times = [], [], []
    for N in orders:
        times = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            MilneSolver(N)
            times.append(time.perf_counter() - t0)
        Ns.append(N)
        mean_times.append(float(np.mean(times)))
        std_times.append(float(np.std(times)))

    return {
        "N": np.array(Ns),
        "mean_time_s": np.array(mean_times),
        "std_time_s": np.array(std_times),
    }


def cli():
    parser = argparse.ArgumentParser(description="Milne PN convergence benchmark")
    parser.add_argument("--nmax", type=int, default=79)
    parser.add_argument("--nmin", type=int, default=1)
    args = parser.parse_args()

    result = convergence_benchmark(n_max=args.nmax, n_min=args.nmin)
    header = "N".rjust(5) + "z0".rjust(16) + "error".rjust(16) + "#modes".rjust(10)
    print(header)
    print("-" * 47)
    for N, z0, err, modes in zip(result["N"], result["z0"], result["error"], result["num_discrete_modes"]):
        print(f"{N:5d}{z0:16.8f}{err:+16.8f}{modes:10d}")


if __name__ == "__main__":
    cli()
