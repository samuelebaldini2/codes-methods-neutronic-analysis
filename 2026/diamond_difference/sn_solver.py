import numpy as np


def gauss_legendre_directions(N):
    if N % 2 != 0:
        raise ValueError("N (number of discrete ordinates) must be even.")
    mu, w = np.polynomial.legendre.leggauss(N)
    return mu, w


def _mirror_pairing(mu):
    N = len(mu)
    pair = np.empty(N, dtype=int)
    for m in range(N):
        pair[m] = np.argmin(np.abs(mu + mu[m]))
    return pair


class Mesh:
    def __init__(self, d, I):
        self.d = d
        self.I = I
        self.dx = np.full(I, d / I)
        self.edges = np.linspace(0.0, d, I + 1)
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])


def make_uniform_source(mesh, sext):
    return np.full(mesh.I, sext)


def make_peaked_source(mesh, sext, width):
    d = mesh.d
    I = mesh.I
    dx = mesh.dx[0]
    # from physical width (e.g. 20cm) to cells: min 1 cell, max whole domain
    n_cells = max(1, int(round(width / dx)))
    n_cells = min(n_cells, I)
    i_mid = I // 2
    # lower bound
    i_lo = max(0, i_mid - n_cells // 2)
    # upper bound
    i_hi = min(I, i_lo + n_cells)
    i_lo = i_hi - n_cells  # re-adjust in case of clipping at hi
    in_region = np.zeros(I, dtype=bool)
    in_region[i_lo:i_hi] = True  # activate source region
    actual_width = n_cells * dx  # in case of clipping this is < the requested width
    Q = np.zeros(I)
    Q[in_region] = sext * d / actual_width
    return Q, actual_width, n_cells


class SNSlabSolver:
    def __init__(self, mesh, sigma_t, c, Q_ext, N=4,
                 bc_left="void", bc_right="void"):
        """
        mesh      : Mesh object
        sigma_t   : scalar or array (I,) total cross section [cm^-1]
        c         : scalar or array (I,) scattering ratio Sigma_s/Sigma_t
        Q_ext     : array (I,) external isotropic source [cm^-3 s^-1]
        N         : quadrature order (even)
        bc_left, bc_right : 'void' or 'reflect'
        """
        self.mesh = mesh
        I = mesh.I
        self.sigma_t = np.full(I, sigma_t) if np.isscalar(sigma_t) else np.asarray(sigma_t, dtype=float)
        self.c = np.full(I, c) if np.isscalar(c) else np.asarray(c, dtype=float)
        self.sigma_s = self.c * self.sigma_t
        self.Q_ext = np.asarray(Q_ext, dtype=float)
        assert bc_left in ("void", "reflect")
        assert bc_right in ("void", "reflect")
        self.bc_left = bc_left
        self.bc_right = bc_right

        self.mu, self.w = gauss_legendre_directions(N)
        self.pair = _mirror_pairing(self.mu)
        self.pos = np.where(self.mu > 0)[0]
        self.neg = np.where(self.mu < 0)[0]

    def _sweep(self, phi_old, psi_left_exit_prev, psi_right_exit_prev):
        mesh = self.mesh
        I = mesh.I
        dx = mesh.dx
        sigma_t = self.sigma_t
        Q_iso = 0.5 * (self.sigma_s * phi_old + self.Q_ext)

        psi_cell = np.zeros((len(self.mu), I))
        psi_left_exit = np.zeros(len(self.mu))   # exiting flux at x=0, per direction (meaningful for mu<0)
        psi_right_exit = np.zeros(len(self.mu))  # exiting flux at x=d, per direction (meaningful for mu>0)

        # positive directions: sweep left -> right
        for m in self.pos:
            mu_m = self.mu[m]
            if self.bc_left == "void":
                psi_edge = 0.0
            else:  # reflect: use mirrored direction's previous exiting flux at left edge
                psi_edge = psi_left_exit_prev[self.pair[m]]
            coef = 2.0 * mu_m / dx
            for i in range(I):
                denom = coef[i] + sigma_t[i]
                psi_i = (Q_iso[i] + coef[i] * psi_edge) / denom
                psi_cell[m, i] = psi_i
                psi_edge = 2.0 * psi_i - psi_edge
            psi_right_exit[m] = psi_edge

        # negative directions: sweep right -> left
        for m in self.neg:
            mu_abs = -self.mu[m]
            if self.bc_right == "void":
                psi_edge = 0.0
            else:  # reflect
                psi_edge = psi_right_exit_prev[self.pair[m]]
            coef = 2.0 * mu_abs / dx
            for i in range(I - 1, -1, -1):
                denom = coef[i] + sigma_t[i]
                psi_i = (Q_iso[i] + coef[i] * psi_edge) / denom
                psi_cell[m, i] = psi_i
                psi_edge = 2.0 * psi_i - psi_edge
            psi_left_exit[m] = psi_edge

        phi_new = self.w @ psi_cell
        return phi_new, psi_cell, psi_left_exit, psi_right_exit

    def solve(self, tol=1e-8, max_iter=20000, verbose=False):
        I = self.mesh.I
        phi_old = np.zeros(I)
        psi_left_exit_prev = np.zeros(len(self.mu))
        psi_right_exit_prev = np.zeros(len(self.mu))

        for it in range(1, max_iter + 1):
            phi_new, psi_cell, psi_left_exit, psi_right_exit = self._sweep(
                phi_old, psi_left_exit_prev, psi_right_exit_prev)

            denom = np.max(np.abs(phi_new))
            denom = denom if denom > 0 else 1.0
            err = np.max(np.abs(phi_new - phi_old)) / denom

            phi_old = phi_new
            psi_left_exit_prev = psi_left_exit
            psi_right_exit_prev = psi_right_exit

            if verbose and (it % 50 == 0 or it == 1):
                print(f"  iter {it:5d}  rel.err = {err:.3e}")

            if err < tol:
                break
        else:
            print(f"WARNING: did not converge in {max_iter} iterations "
                  f"(last err = {err:.3e})")

        self.phi = phi_old
        self.psi_cell = psi_cell
        self.n_iter = it
        self.final_err = err
        return phi_old, it
