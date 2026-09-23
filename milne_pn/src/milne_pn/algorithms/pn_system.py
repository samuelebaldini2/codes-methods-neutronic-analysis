"""
pn_system.py
------------
Builds the matrices needed for the PN approximation to the (purely
scattering, isotropic) Milne problem:

    A * Phi(x) + D * Phi(x) = 0

with A tridiagonal / zero-diagonal:
    A[l, l-1] = l / (2l+1)
    A[l, l+1] = (l+1) / (2l+1)
and D = diag(0, 1, 1, ..., 1)  (size N+1, N odd).

H[k, l] = integral_0^1 P_k(mu) P_l(mu) dmu  is the half-range Legendre
overlap matrix needed for the Marshak boundary conditions; it is computed
exactly with numpy Gauss-Legendre quadrature.
"""

import numpy as np

from ..utils.quadrature import gauss_legendre_half_range, legendre_matrix


class PNSystem:
    def __init__(self, N, quad_multiplier=4):
        if N < 1 or N % 2 == 0:
            raise ValueError("N must be odd and >= 1")
        self.N = N
        self.size = N + 1
        self.A, self.D = self._build_AD()
        self.H = self._build_H(quad_multiplier)

    def _build_AD(self):
        size = self.size
        A = np.zeros((size, size))
        D = np.diag([0.0] + [1.0] * self.N)
        for l in range(size):
            if l - 1 >= 0:
                A[l, l - 1] = l / (2 * l + 1)
            if l + 1 < size:
                A[l, l + 1] = (l + 1) / (2 * l + 1)
        return A, D

    def _build_H(self, quad_multiplier):
        size = self.size
        n_quad = quad_multiplier * size
        mu, wq = gauss_legendre_half_range(n_quad)
        P = legendre_matrix(self.N, mu)
        H = (P * wq[None, :]) @ P.T
        return H
