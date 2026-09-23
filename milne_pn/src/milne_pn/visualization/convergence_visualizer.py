"""
convergence_visualizer.py
--------------------------
z0(N) vs N, and the convergence error on a log scale.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_convergence(N, z0, z0_exact, out_path):
    N = np.asarray(N)
    z0 = np.asarray(z0)
    err = np.abs(z0 - z0_exact)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.axhline(z0_exact, color="crimson", ls="--", lw=1.3, label=f"Case exact: {z0_exact}")
    ax.plot(N, z0, "o-", color="steelblue", ms=4, label=r"$P_N$ result")
    ax.set_xlabel("PN order  N")
    ax.set_ylabel(r"extrapolation distance $z_0$  [mfp]")
    ax.set_title(r"Convergence of $z_0$ with $P_N$ order")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.semilogy(N, np.maximum(err, 1e-16), "o-", color="darkorange", ms=4)
    ax.set_xlabel("PN order  N")
    ax.set_ylabel(r"$|z_0^{(N)} - z_0^{exact}|$")
    ax.set_title("Convergence error (log scale)")
    ax.grid(alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
