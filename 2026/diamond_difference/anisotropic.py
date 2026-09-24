import numpy as np

from sn_solver import SNSlabSolver  # reused as-is, not modified

_G_MAX = 1.0 / 3.0  # |g| <= 1/3 keeps the phase function non-negative


def phase_function(mu0, g):
    return 0.5 * (1.0 + 3.0 * g * mu0)


def scattering_moments(sigma_t, c, g):
    sigma_s0 = c * sigma_t
    sigma_s1 = g * sigma_s0
    return sigma_s0, sigma_s1


def source_angular_weight(mu, shape):
    if shape == "isotropic":
        return 0.5
    elif shape == "forward_peaked":
        return 0.5 * (1.0 + mu)
    raise ValueError(f"Unknown source angular shape: '{shape}' "
                      f"(expected 'isotropic' or 'forward_peaked')")


def direction_source(mu, phi0, phi1, sigma_s0, sigma_s1, Q_ext, source_shape):
    return (0.5 * sigma_s0 * phi0
            + 1.5 * sigma_s1 * mu * phi1
            + Q_ext * source_angular_weight(mu, source_shape))


class AnisotropicSNSlabSolver(SNSlabSolver):
    def __init__(self, mesh, sigma_t, c, Q_ext, N=4,
                 bc_left="void", bc_right="void",
                 g=0.0, source_shape="isotropic"):
        super().__init__(mesh, sigma_t, c, Q_ext, N=N,
                          bc_left=bc_left, bc_right=bc_right)
        if source_shape not in ("isotropic", "forward_peaked"):
            raise ValueError(f"Unknown source_shape: '{source_shape}'")
        g_arr = np.full(mesh.I, g) if np.isscalar(g) else np.asarray(g, dtype=float)
        if np.any(np.abs(g_arr) > _G_MAX):
            raise ValueError(f"anisotropy factor g must satisfy |g| <= {_G_MAX:.4f} "
                              f"for the phase function to stay non-negative")
        # self.sigma_s (= Sigma_s,0 = c*Sigma_t) already set by SNSlabSolver.__init__
        self.g = g_arr
        _, self.sigma_s1 = scattering_moments(self.sigma_t, self.c, self.g)
        self.source_shape = source_shape

    #sweep but the source is direction dependent (rebuilt once per direction)
    def _sweep(self, phi0_old, phi1_old, psi_left_exit_prev, psi_right_exit_prev):
        mesh = self.mesh
        I = mesh.I
        dx = mesh.dx
        sigma_t = self.sigma_t
        psi_cell = np.zeros((len(self.mu), I))
        psi_left_exit = np.zeros(len(self.mu))
        psi_right_exit = np.zeros(len(self.mu))

        for m in self.pos:
            mu_m = self.mu[m]
            Q_dir = direction_source(mu_m, phi0_old, phi1_old, self.sigma_s, self.sigma_s1,self.Q_ext, self.source_shape)
            psi_edge = 0.0 if self.bc_left == "void" else psi_left_exit_prev[self.pair[m]]
            coef = 2.0 * mu_m / dx
            for i in range(I):
                denom = coef[i] + sigma_t[i]
                psi_i = (Q_dir[i] + coef[i] * psi_edge) / denom
                psi_cell[m, i] = psi_i
                psi_edge = 2.0 * psi_i - psi_edge
            psi_right_exit[m] = psi_edge

        for m in self.neg:
            mu_m = self.mu[m]  # negative; used as-is (signed) in direction_source
            mu_abs = -mu_m
            Q_dir = direction_source(mu_m, phi0_old, phi1_old,
                                      self.sigma_s, self.sigma_s1,
                                      self.Q_ext, self.source_shape)
            psi_edge = 0.0 if self.bc_right == "void" else psi_right_exit_prev[self.pair[m]]
            coef = 2.0 * mu_abs / dx
            for i in range(I - 1, -1, -1):
                denom = coef[i] + sigma_t[i]
                psi_i = (Q_dir[i] + coef[i] * psi_edge) / denom
                psi_cell[m, i] = psi_i
                psi_edge = 2.0 * psi_i - psi_edge
            psi_left_exit[m] = psi_edge

        phi0_new = self.w @ psi_cell
        phi1_new = (self.w * self.mu) @ psi_cell
        return phi0_new, phi1_new, psi_cell, psi_left_exit, psi_right_exit

    def solve(self, tol=1e-8, max_iter=20000, verbose=False):
        I = self.mesh.I
        phi0_old = np.zeros(I)
        phi1_old = np.zeros(I)
        psi_left_exit_prev = np.zeros(len(self.mu))
        psi_right_exit_prev = np.zeros(len(self.mu))

        for it in range(1, max_iter + 1):
            phi0_new, phi1_new, psi_cell, psi_left_exit, psi_right_exit = self._sweep(
                phi0_old, phi1_old, psi_left_exit_prev, psi_right_exit_prev)

            denom = np.max(np.abs(phi0_new))
            denom = denom if denom > 0 else 1.0
            err = np.max(np.abs(phi0_new - phi0_old)) / denom

            phi0_old, phi1_old = phi0_new, phi1_new
            psi_left_exit_prev = psi_left_exit
            psi_right_exit_prev = psi_right_exit

            if verbose and (it % 50 == 0 or it == 1):
                print(f"  iter {it:5d}  rel.err = {err:.3e}")

            if err < tol:
                break
        else:
            print(f"WARNING: did not converge in {max_iter} iterations "
                  f"(last err = {err:.3e})")

        self.phi = phi0_old   
        self.phi0 = phi0_old  
        self.phi1 = phi1_old  
        self.psi_cell = psi_cell
        self.n_iter = it
        self.final_err = err
        return phi0_old, it
