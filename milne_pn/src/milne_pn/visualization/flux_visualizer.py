"""
flux_visualizer.py
-------------------
phi0(x), its linear asymptote x + z0, and the extrapolated zero-crossing.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_flux_profile(x, phi0, z0, out_path, N_label=None):
    x = np.asarray(x)
    phi0 = np.asarray(phi0)
    phi0_asym = x + z0

    fig, ax = plt.subplots(figsize=(6.5, 5))
    label_full = rf"$\phi_0(x)$ (full $P_{{{N_label}}}$ solution)" if N_label else r"$\phi_0(x)$ (full $P_N$ solution)"
    ax.plot(x, phi0, color="steelblue", lw=2, label=label_full)
    ax.plot(x, phi0_asym, color="crimson", ls="--", lw=1.5, label=r"asymptotic term: $x + z_0$")

    xline = np.linspace(-z0, 0, 50)
    ax.plot(xline, xline + z0, color="crimson", ls=":", lw=1.5)
    ax.axvline(0, color="gray", lw=0.8)
    ax.axhline(0, color="gray", lw=0.8)
    ax.scatter([-z0], [0], color="crimson", zorder=5)
    ax.annotate(rf"$-z_0={-z0:.4f}$", (-z0, 0), textcoords="offset points", xytext=(8, -14))

    ax.set_xlim(-1, max(6.0, float(x.max())))
    ax.set_xlabel("x  [mean free paths]")
    ax.set_ylabel(r"$\phi_0(x)$")
    title = rf"Milne problem scalar flux, $P_{{{N_label}}}$ approximation" if N_label else "Milne problem scalar flux"
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
