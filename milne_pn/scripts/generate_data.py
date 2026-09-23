#!/usr/bin/env python3
"""
scripts/generate_data.py
--------------------------
Thin CLI wrapper: generates all the CSV data products used by the
visualization layer and the notebooks -- convergence.csv, flux_profile.csv,
and gibbs_overshoot.csv -- into data/output/.

Usage:
    python scripts/generate_data.py
    python scripts/generate_data.py --nmax 99 --nplot 41
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from milne_pn.algorithms.milne_solver import MilneSolver
from milne_pn.algorithms.runge_analysis import gibbs_convergence_table
from milne_pn.benchmarks.benchmark import convergence_benchmark
from milne_pn.utils.config import load_config
from milne_pn.utils.io_utils import write_csv


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, "..")
    cfg = load_config(os.path.join(root, "config", "config.yaml"))

    parser = argparse.ArgumentParser(description="Generate Milne PN data products")
    parser.add_argument("--nmax", type=int, default=cfg["solver"]["n_max"])
    parser.add_argument("--nplot", type=int, default=cfg["solver"]["n_plot"])
    parser.add_argument("--out-dir", type=str, default=os.path.join(root, cfg["paths"]["data_output"]))
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    z0_exact = cfg["reference"]["z0_exact"]

    conv = convergence_benchmark(n_max=args.nmax, z0_exact=z0_exact)
    write_csv(os.path.join(args.out_dir, "convergence.csv"), {
        "N": conv["N"], "z0": conv["z0"], "error": conv["error"],
        "num_discrete_modes": conv["num_discrete_modes"],
        "max_imag_residual": conv["max_imag_residual"],
    })
    print(f"Wrote {os.path.join(args.out_dir, 'convergence.csv')}")

    solver = MilneSolver(args.nplot)
    x = np.linspace(0.0, 6.0, 401)
    phi0 = solver.phi0(x)
    write_csv(os.path.join(args.out_dir, "flux_profile.csv"), {
        "x": x, "phi0": phi0, "phi0_asymptotic": x + solver.z0,
    })
    print(f"Wrote {os.path.join(args.out_dir, 'flux_profile.csv')} "
          f"(P{args.nplot}, z0={solver.z0:.6f})")

    orders = cfg["runge_analysis"]["orders"]
    rows = gibbs_convergence_table(orders, n_mu=cfg["runge_analysis"]["n_mu"])
    write_csv(os.path.join(args.out_dir, "gibbs_overshoot.csv"), {
        "N": np.array([r["N"] for r in rows]),
        "max_overshoot_abs": np.array([r["max_overshoot_abs"] for r in rows]),
        "outgoing_scale": np.array([r["outgoing_scale"] for r in rows]),
        "relative_overshoot": np.array([r["relative_overshoot"] for r in rows]),
    })
    print(f"Wrote {os.path.join(args.out_dir, 'gibbs_overshoot.csv')}")


if __name__ == "__main__":
    main()
