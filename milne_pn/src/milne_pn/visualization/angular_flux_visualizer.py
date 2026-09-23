"""
angular_flux_visualizer.py
----------------------------
Visualizes the Gibbs ringing in the PN-reconstructed boundary angular flux
psi(0, mu), and the classical Runge phenomenon as a side-by-side analogue.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ..algorithms.milne_solver import MilneSolver
from ..algorithms.runge_analysis import (
    reconstruct_boundary_flux,
    gibbs_overshoot,
    runge_interpolation_demo,
)


def plot_gibbs_ringing(orders, out_path, n_mu=2000):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(orders)))
    for N, color in zip(orders, colors):
        solver = MilneSolver(N)
        mu, psi0 = reconstruct_boundary_flux(solver, n_mu=n_mu)
        ax.plot(mu, psi0, color=color, lw=1.4, label=rf"$P_{{{N}}}$")
    ax.axvline(0, color="black", lw=0.8)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvspan(0, 1, color="red", alpha=0.05)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\psi(0,\mu)$  (truncated PN reconstruction)")
    ax.set_title("Gibbs ringing at the angular discontinuity ($\\mu=0$)\n"
                  "(shaded: should be exactly 0 -- vacuum b.c.)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    Ns = np.array(orders)
    rel = np.array([gibbs_overshoot(MilneSolver(N), n_mu=n_mu)["relative_overshoot"] for N in orders])
    ax.plot(Ns, rel * 100, "o-", color="firebrick")
    ax.set_xlabel("PN order  N")
    ax.set_ylabel("max overshoot in $\\mu>0$,\nrelative to outgoing flux scale [%]")
    ax.set_title("Overshoot amplitude does NOT vanish with N\n(classic Gibbs plateau)")
    ax.grid(alpha=0.3)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_runge_phenomenon(degrees, out_path, n_eval=1000):
    result = runge_interpolation_demo(degrees=degrees, n_eval=n_eval)
    x = result["x_eval"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    ax = axes[0]
    ax.plot(x, result["y_true"], "k-", lw=2, label=r"$f(x)=\dfrac{1}{1+25x^2}$")
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(degrees)))
    for d, color in zip(degrees, colors):
        ax.plot(x, result["interpolants"][d], color=color, lw=1.2, label=f"degree {d}")
    ax.set_ylim(-1, 2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Classical Runge phenomenon\n(equispaced-node polynomial interpolation)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    degs = sorted(result["max_abs_error"].keys())
    errs = [result["max_abs_error"][d] for d in degs]
    ax.semilogy(degs, errs, "o-", color="firebrick")
    ax.set_xlabel("interpolation degree")
    ax.set_ylabel("max abs. interpolation error")
    ax.set_title("Error GROWS with degree\n(pathology of global polynomials, not noise)")
    ax.grid(alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
