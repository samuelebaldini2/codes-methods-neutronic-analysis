"""
milne_solver.py
----------------
Assembles and solves the PN Milne problem for a given odd order N.

  1. Forms M = A^{-1} D and finds its eigenvalues/eigenvectors (numpy
     built-ins: numpy.linalg.solve, numpy.linalg.eig) -- the exponential
     modes Phi(x) = v * exp(-omega x).
  2. Keeps only the physically bounded (Re(omega) > 0, decaying as
     x -> +infinity) discrete modes, plus the exact analytic asymptotic
     (linear-in-x diffusion) Jordan-chain solution associated with the
     defective double omega = 0 eigenvalue.
  3. Applies Marshak boundary conditions at x = 0 to fix the mode
     amplitudes, including P0 -- which IS the extrapolation distance z0
     by construction (numpy.linalg.solve; no numerical fit involved).
"""

import numpy as np

from .pn_system import PNSystem


class MilneSolver:
    def __init__(self, N, tol=1e-9, quad_multiplier=4):
        self.N = N
        self.tol = tol
        self.pn = PNSystem(N, quad_multiplier=quad_multiplier)
        self._solve()

    def _solve(self):
        pn = self.pn
        size = pn.size

        M = np.linalg.solve(pn.A, pn.D)
        w, V = np.linalg.eig(M)

        self.max_imag_residual = float(np.max(np.abs(w.imag)))

        nonzero = np.abs(w) > self.tol
        w_nz, V_nz = w[nonzero], V[:, nonzero]

        decay = w_nz.real > self.tol
        w_dec, V_dec = w_nz[decay].real, V_nz[:, decay].real
        n_modes = w_dec.size

        e0 = np.zeros(size); e0[0] = 1.0
        w_vec = np.zeros(size)
        if size > 1:
            w_vec[1] = -1.0 / 3.0

        H = pn.H
        odd_ls = list(range(1, self.N + 1, 2))
        n_unknown = 1 + n_modes
        assert len(odd_ls) == n_unknown, (len(odd_ls), n_unknown)

        coeffs_l = np.array([(2 * l + 1) / 2 for l in range(size)])

        Mtx = np.zeros((len(odd_ls), n_unknown))
        rhs = np.zeros(len(odd_ls))

        for irow, k in enumerate(odd_ls):
            Hk = H[k, :]
            Mtx[irow, 0] = np.sum(coeffs_l * Hk * e0)
            for i in range(n_modes):
                Mtx[irow, 1 + i] = np.sum(coeffs_l * Hk * V_dec[:, i])
            rhs[irow] = -np.sum(coeffs_l * Hk * w_vec)

        sol = np.linalg.solve(Mtx, rhs)
        self.P0 = sol[0]
        self.c = sol[1:]
        self.omega = w_dec
        self.V = V_dec

        self.z0 = self.P0

    def phi0(self, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        val = self.P0 + x
        for i in range(self.omega.size):
            val = val + self.c[i] * self.V[0, i] * np.exp(-self.omega[i] * x)
        return val

    def phi_l(self, l, x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        if l == 0:
            return self.phi0(x)
        if l == 1:
            val = np.full_like(x, -1.0 / 3.0)
        else:
            val = np.zeros_like(x)
        for i in range(self.omega.size):
            val = val + self.c[i] * self.V[l, i] * np.exp(-self.omega[i] * x)
        return val

    def angular_flux(self, x, mu):
        from ..utils.quadrature import legendre_matrix
        x = np.atleast_1d(np.asarray(x, dtype=float))
        mu = np.atleast_1d(np.asarray(mu, dtype=float))
        P = legendre_matrix(self.N, mu)
        psi = np.zeros((x.size, mu.size))
        for l in range(self.N + 1):
            coeff = (2 * l + 1) / 2.0
            phil = self.phi_l(l, x)
            psi += coeff * np.outer(phil, P[l])
        return psi

    @property
    def num_discrete_modes(self):
        return int(self.omega.size)
