#!/usr/bin/env python3
"""
scripts/benchmark.py
---------------------
Thin CLI wrapper: run the convergence benchmark using config/config.yaml
and write results to data/output/convergence.csv (plus a timing table if
--timing is passed).

Usage:
    python scripts/benchmark.py
    python scripts/benchmark.py --nmax 99
    python scripts/benchmark.py --timing
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from milne_pn.benchmarks.benchmark import convergence_benchmark, timing_benchmark
from milne_pn.utils.config import load_config
from milne_pn.utils.io_utils import write_csv


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.join(here, "..")
    cfg = load_config(os.path.join(root, "config", "config.yaml"))

    parser = argparse.ArgumentParser(description="Milne PN convergence/timing benchmark")
    parser.add_argument("--nmax", type=int, default=cfg["solver"]["n_max"])
    parser.add_argument("--nmin", type=int, default=1)
    parser.add_argument("--timing", action="store_true", help="also run the timing benchmark")
    parser.add_argument("--out-dir", type=str, default=os.path.join(root, cfg["paths"]["data_output"]))
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    result = convergence_benchmark(n_max=args.nmax, n_min=args.nmin, z0_exact=cfg["reference"]["z0_exact"])
    header = "N".rjust(5) + "z0".rjust(16) + "error".rjust(16) + "#modes".rjust(10)
    print(header)
    print("-" * 47)
    for N, z0, err, modes in zip(result["N"], result["z0"], result["error"], result["num_discrete_modes"]):
        print(f"{N:5d}{z0:16.8f}{err:+16.8f}{modes:10d}")

    conv_csv = os.path.join(args.out_dir, "convergence.csv")
    write_csv(conv_csv, {
        "N": result["N"], "z0": result["z0"], "error": result["error"],
        "num_discrete_modes": result["num_discrete_modes"],
        "max_imag_residual": result["max_imag_residual"],
    })
    print(f"\nWrote {conv_csv}")

    if args.timing:
        orders = list(range(5, args.nmax + 1, 10))
        timing = timing_benchmark(orders)
        timing_csv = os.path.join(args.out_dir, "timing.csv")
        write_csv(timing_csv, timing)
        print(f"Wrote {timing_csv}")


if __name__ == "__main__":
    main()
