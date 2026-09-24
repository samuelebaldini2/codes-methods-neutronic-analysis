import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  


def _save(fig, out_path):
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_polar_polynomial_exactness(n_mu, rows, out_path):
    k = np.array([r[0] for r in rows])
    err = np.array([r[3] for r in rows])
    within = np.array([r[4] for r in rows])
    err_plot = np.maximum(err, 1e-17)  

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(k[within], err_plot[within], "o-", color="tab:blue",
                label="within guaranteed range")
    ax.semilogy(k[~within], err_plot[~within], "o-", color="tab:red",
                label="beyond guaranteed range")
    ax.axvline(2 * n_mu - 1, color="gray", linestyle="--", linewidth=1,
               label=f"degree limit = 2*n_mu-1 = {2 * n_mu - 1}")
    ax.set_xlabel("degree k of mu^k")
    ax.set_ylabel("absolute error |numeric - exact|")
    ax.set_title(f"Polynomial exactness in mu (Gauss-Legendre, n_mu={n_mu})")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    return _save(fig, out_path)


def plot_spherical_harmonic_orthonormality(n_mu, n_phi, azimuthal, n_max,
                                            err_records, out_path):
    errs = np.array([r[4] for r in err_records])
    diag = np.array([r[5] for r in err_records])
    err_plot = np.maximum(errs, 1e-17)
    order = np.argsort(-err_plot)  
    idx = np.arange(len(err_records))

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogy(idx[~diag[order]], err_plot[order][~diag[order]], "o",
                color="tab:orange", markersize=3, label="off-diagonal pairs (n1,m1)!=(n2,m2)")
    ax.semilogy(idx[diag[order]], err_plot[order][diag[order]], "o",
                color="tab:blue", markersize=3, label="diagonal pairs (n1,m1)=(n2,m2)")
    ax.axhline(1e-10, color="gray", linestyle="--", linewidth=1,
               label="check threshold 1e-10")
    ax.set_xlabel(f"pairs (n1,m1,n2,m2), n1,n2<={n_max} -- sorted by decreasing error")
    ax.set_ylabel("absolute error")
    ax.set_title(f"Spherical harmonic orthonormality\n"
                 f"n_mu={n_mu}, n_phi={n_phi}, azimuthal='{azimuthal}'")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    return _save(fig, out_path)


def plot_orthonormality_beyond_resolution(n_mu, n_phi, curve, n_too_high, err, out_path):
    n = np.array([c[0] for c in curve])
    err_n = np.maximum(np.array([c[1] for c in curve]), 1e-17)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(n, err_n, "o-", color="tab:blue")
    ax.axvline(n_mu, color="gray", linestyle="--", linewidth=1,
               label=f"n_mu={n_mu} (expected resolution limit)")
    ax.axhline(1e-6, color="tab:red", linestyle=":", linewidth=1,
               label="'appreciable error' threshold 1e-6")
    ax.scatter([n_too_high], [max(err, 1e-17)], color="tab:red", zorder=5,
               label=f"n tested in the check = {n_too_high}")
    ax.set_xlabel("n (degree of Y_n^n tested against itself)")
    ax.set_ylabel("absolute error")
    ax.set_title(f"Quadrature resolution limit (n_mu={n_mu}, n_phi={n_phi})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    return _save(fig, out_path)


def plot_quadrature_nodes_3d(n_mu, n_phi, azimuthal, x, y, z, color_values, color_label, out_path):
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="lightgray", alpha=0.15, linewidth=0,
                     rstride=1, cstride=1, shade=False)
    if np.any(color_values < 0):
        vmax = np.max(np.abs(color_values))
        cmap, vmin, vmax_plot = "RdBu_r", -vmax, vmax
    else:
        cmap, vmin, vmax_plot = "viridis", None, None

    sc = ax.scatter(x, y, z, c=color_values, cmap=cmap, vmin=vmin, vmax=vmax_plot,
                     s=25, depthshade=True)
    fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, label=color_label)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z  (= mu)")
    ax.set_title(f"Product quadrature nodes on the unit sphere\n"
                 f"n_mu={n_mu}, n_phi={n_phi}, azimuthal='{azimuthal}' "
                 f"({n_mu * n_phi} nodes)")
    ax.set_box_aspect([1, 1, 1])
    return _save(fig, out_path)


def plot_gauss_chebyshev_trapezoidal_equivalence(n_phi, data, out_path):
    phi_trap = data["phi_trap"]
    phi_cheb = data["phi_cheb"]
    k_rows = data["k_rows"]
    k = np.array([r[0] for r in k_rows])
    err_trap = np.maximum(np.array([r[3] for r in k_rows]), 1e-17)
    err_cheb = np.maximum(np.array([r[4] for r in k_rows]), 1e-17)
    within = np.array([r[5] for r in k_rows])

    fig = plt.figure(figsize=(11, 4.8))
    ax_nodes = fig.add_subplot(1, 2, 1, projection="polar")
    ax_err = fig.add_subplot(1, 2, 2)

    ax_nodes.scatter(phi_trap, np.ones_like(phi_trap), color="tab:blue", s=40,
                      label="trapezoidal nodes")
    ax_nodes.scatter(phi_cheb, np.ones_like(phi_cheb), color="tab:orange", s=40,
                      marker="x", label="Gauss-Chebyshev nodes")
    ax_nodes.set_yticklabels([])
    ax_nodes.set_title(f"Azimuthal nodes (n_phi={n_phi})", pad=20)
    ax_nodes.legend(loc="upper right", bbox_to_anchor=(1.3, 1.15), fontsize=8)

    ax_err.semilogy(k[within], err_trap[within], "o-", color="tab:blue",
                     label="trapezoidal")
    ax_err.semilogy(k[within], err_cheb[within], "x--", color="tab:orange",
                     label="Gauss-Chebyshev")
    ax_err.semilogy(k[~within], err_trap[~within], "o-", color="tab:blue", alpha=0.4)
    ax_err.semilogy(k[~within], err_cheb[~within], "x--", color="tab:orange", alpha=0.4)
    ax_err.axvline(n_phi - 1, color="gray", linestyle="--", linewidth=1,
                    label=f"degree limit = n_phi-1 = {n_phi - 1}")
    ax_err.set_xlabel("k in e^(i k phi)")
    ax_err.set_ylabel("absolute error")
    ax_err.set_title("Error vs degree k\n(solid = within range, faded = aliasing expected)")
    ax_err.legend(fontsize=8)
    ax_err.grid(True, which="both", alpha=0.3)

    fig.suptitle("Gauss-Chebyshev <-> trapezoidal equivalence")
    return _save(fig, out_path)