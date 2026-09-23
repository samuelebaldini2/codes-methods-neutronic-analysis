"""
runge_analysis.py
------------------
The exact boundary angular flux for the Milne problem has a jump
discontinuity at mu = 0: zero for mu > 0 (vacuum), nonzero for mu < 0
(outgoing). The PN method reconstructs psi(0, mu) as a finite global
Legendre series, which is exactly the setting of the GIBBS PHENOMENON:
the reconstruction rings near the discontinuity with an amplitude that
does NOT shrink as N grows, only the ringing width does. This is the
angular-domain sibling of the classical RUNGE PHENOMENON (global
polynomial interpolation ringing near domain edges as degree grows).
Both explain why PN convergence of z0 is algebraic, not exponential.
"""

import numpy as np

from ..algorithms.milne_solver import MilneSolver


def reconstruct_boundary_flux(solver, n_mu=2000):
    mu = np.linspace(-1.0, 1.0, n_mu)
    psi0 = solver.angular_flux(0.0, mu)[0]
    return mu, psi0


def gibbs_overshoot(solver, n_mu=2000):
    mu, psi0 = reconstruct_boundary_flux(solver, n_mu)
    incoming = mu > 0
    outgoing = mu < 0

    outgoing_scale = float(np.max(np.abs(psi0[outgoing]))) if np.any(outgoing) else float("nan")
    max_overshoot = float(np.max(np.abs(psi0[incoming]))) if np.any(incoming) else 0.0
    rel_overshoot = max_overshoot / outgoing_scale if outgoing_scale else float("nan")

    return {
        "N": solver.N,
        "max_overshoot_abs": max_overshoot,
        "outgoing_scale": outgoing_scale,
        "relative_overshoot": rel_overshoot,
    }


def gibbs_convergence_table(orders, n_mu=2000):
    rows = []
    for N in orders:
        solver = MilneSolver(N)
        rows.append(gibbs_overshoot(solver, n_mu=n_mu))
    return rows


def runge_function(x):
    x = np.asarray(x, dtype=float)
    return 1.0 / (1.0 + 25.0 * x ** 2)


def runge_interpolation_demo(degrees=(5, 10, 15, 20), n_eval=1000):
    x_eval = np.linspace(-1.0, 1.0, n_eval)
    y_true = runge_function(x_eval)

    interpolants = {}
    max_abs_error = {}
    for d in degrees:
        x_nodes = np.linspace(-1.0, 1.0, d + 1)
        y_nodes = runge_function(x_nodes)
        coeffs = np.polyfit(x_nodes, y_nodes, d)
        y_interp = np.polyval(coeffs, x_eval)
        interpolants[d] = y_interp
        max_abs_error[d] = float(np.max(np.abs(y_interp - y_true)))

    return {
        "x_eval": x_eval,
        "y_true": y_true,
        "interpolants": interpolants,
        "max_abs_error": max_abs_error,
    }
