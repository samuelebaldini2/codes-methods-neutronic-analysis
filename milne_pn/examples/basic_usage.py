#!/usr/bin/env python3
"""
examples/basic_usage.py
-------------------------
The smallest complete example: solve the Milne problem at a single PN
order and read off the extrapolation distance z0.
"""

from milne_pn import MilneSolver


def main():
    N = 21
    solver = MilneSolver(N)

    print(f"PN order N = {N}")
    print(f"z0 = {solver.z0:.6f} mean free paths  (Case exact: 0.710446)")
    print(f"number of discrete decaying modes kept: {solver.num_discrete_modes}")

    print("\nScalar flux phi0(x) at a few depths:")
    for x in (0.0, 0.5, 1.0, 2.0, 5.0):
        print(f"  x = {x:5.2f}   phi0(x) = {solver.phi0(x)[0]:.6f}")


if __name__ == "__main__":
    main()
