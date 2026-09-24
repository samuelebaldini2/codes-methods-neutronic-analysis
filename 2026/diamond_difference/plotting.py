
import numpy as np
import matplotlib.pyplot as plt

from anisotropic import phase_function, source_angular_weight

_PALETTE = ["tab:blue", "tab:red", "tab:green", "tab:orange", "tab:purple", "tab:brown"]


def plot_anisotropic_functions(out_dir, g, sext, source_shape):
    mu = np.linspace(-1.0, 1.0, 400)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.plot(mu, phase_function(mu, g), lw=2, color="tab:blue", label=f"g={g}")
    ax.axhline(0.5, ls="--", color="gray", lw=1, label="isotropic (g=0)")
    ax.set_xlabel(r"$\mu_0$")
    ax.set_ylabel(r"$p(\mu_0)$")
    ax.set_title("Scattering phase function")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    s_ext = sext * np.array([source_angular_weight(m, source_shape) for m in mu])
    ax.plot(mu, s_ext, lw=2, color="tab:red", label=source_shape)
    ax.axhline(0.5 * sext, ls="--", color="gray", lw=1, label="isotropic")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$s_{ext}(\mu)$")
    ax.set_title("External source angular shape")
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    path = f"{out_dir}/anisotropic_functions.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_bc_comparison(out_dir, curves, title):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for i, (label, mesh, phi) in enumerate(curves):
        ax.plot(mesh.centers, phi, lw=2, color=_PALETTE[i % len(_PALETTE)], label=label)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(r"$\phi(x)$")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = f"{out_dir}/solution.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_solution(out_dir, mesh, phi, nit, title, phi0_analytic=None):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(mesh.centers, phi, lw=2, color="tab:blue", label=f"SN solution ({nit} it.)")
    if phi0_analytic is not None:
        ax.axhline(phi0_analytic, ls="--", color="k",
                    label=f"analytical $\\phi_0={phi0_analytic:.3f}$")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(r"$\phi(x)$")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = f"{out_dir}/solution.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_grid_convergence(out_dir, gc_cell_counts, min_edges, plot_results):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ax = axes[0]
    ax.semilogx(gc_cell_counts, min_edges, "o-", color="tab:purple")
    ax.axhline(0, color="k", lw=1)
    ax.set_xlabel("number of cells I")
    ax.set_ylabel("min. edge angular flux over mesh & directions")
    ax.set_title("Diamond-difference positivity vs mesh refinement")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for I, (mesh, phi) in plot_results.items():
        ax.plot(mesh.centers, phi, label=f"I={I}")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(r"$\phi(x)$")
    ax.set_title("Scalar flux convergence with mesh refinement")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    path = f"{out_dir}/grid_convergence.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_iterations_vs_c(out_dir, curves, tol, title, show_legend):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for i, (label, c_values, iters) in enumerate(curves):
        ax.semilogy(c_values, iters, "o-", ms=3, color=_PALETTE[i % len(_PALETTE)], label=label)
    ax.set_xlabel("scattering ratio c")
    ax.set_ylabel(f"number of source-iterations to reach tol={tol:.0e}")
    ax.set_title(title)
    if show_legend:
        ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    path = f"{out_dir}/iterations_vs_c.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
