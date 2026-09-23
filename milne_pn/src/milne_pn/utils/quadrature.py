"""
quadrature.py
-------------
Numerical-quadrature helpers built on numpy own routines
(numpy.polynomial.legendre), used for the half-range Legendre overlap
integrals that appear in the Marshak boundary conditions, and for
reconstructing/plotting truncated Legendre series.
"""

import numpy as np


def gauss_legendre_half_range(n_quad: int):
    """Gauss-Legendre nodes/weights mapped from [-1, 1] onto the half range
    [0, 1] (used for integral_0^1 f(mu) dmu, as in the Marshak conditions)."""
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_quad)
    mu = 0.5 * (x_gl + 1.0)
    wq = 0.5 * w_gl
    return mu, wq


def legendre_matrix(max_degree: int, mu: np.ndarray) -> np.ndarray:
    """Return P[l, i] = P_l(mu[i]) for l = 0..max_degree, via
    numpy.polynomial.legendre.legval."""
    mu = np.asarray(mu, dtype=float)
    P = np.zeros((max_degree + 1, mu.size))
    for l in range(max_degree + 1):
        c = np.zeros(l + 1)
        c[-1] = 1.0
        P[l] = np.polynomial.legendre.legval(mu, c)
    return P
