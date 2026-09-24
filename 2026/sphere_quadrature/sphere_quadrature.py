import numpy as np


def gauss_legendre_polar(n_mu):
    mu, w_mu = np.polynomial.legendre.leggauss(n_mu)
    theta = np.arccos(mu)
    return mu, w_mu, theta


def trapezoidal_azimuthal(n_phi, phi0=0.0):
    phi = phi0 + 2 * np.pi * np.arange(n_phi) / n_phi
    w_phi = np.full(n_phi, 2 * np.pi / n_phi)
    return phi, w_phi


def gauss_chebyshev_azimuthal(n_phi):
    k = np.arange(n_phi)
    phi = (2 * k + 1) * np.pi / n_phi
    w_phi = np.full(n_phi, 2 * np.pi / n_phi)
    return phi, w_phi


AZIMUTHAL_RULES = {
    "trapezoidal": trapezoidal_azimuthal,
    "gauss_chebyshev": gauss_chebyshev_azimuthal,
}


def sphere_quadrature(n_mu, n_phi, azimuthal="trapezoidal"):
    if azimuthal not in AZIMUTHAL_RULES:
        raise ValueError(f"azimuthal must be one of {list(AZIMUTHAL_RULES)}")

    mu, w_mu, theta = gauss_legendre_polar(n_mu)
    phi, w_phi = AZIMUTHAL_RULES[azimuthal](n_phi)

    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(w_mu, w_phi)

    return {
        "mu": mu, "w_mu": w_mu, "theta": theta,
        "phi": phi, "w_phi": w_phi,
        "THETA": THETA, "PHI": PHI, "W": W,
    }


def integrate_sphere(f, quad):
    return np.sum(f(quad["THETA"], quad["PHI"]) * quad["W"])
