import numpy as np

from sph_harm_compat import sph_harm_y
from sphere_quadrature import (
    gauss_legendre_polar, trapezoidal_azimuthal, gauss_chebyshev_azimuthal,
    sphere_quadrature, integrate_sphere,
)


def analytic_mu_power_integral(k):
    if k % 2 == 1:
        return 0.0
    return 2.0 / (k + 1)


def polar_polynomial_errors(n_mu, k_max=None):
    if k_max is None:
        k_max = 2 * n_mu + 2
    mu, w_mu, _ = gauss_legendre_polar(n_mu)

    rows = []
    max_err_within = 0.0
    for k in range(0, k_max + 1):
        numeric = np.sum(w_mu * mu ** k)
        exact = analytic_mu_power_integral(k)
        err = abs(numeric - exact)
        within_range = k <= 2 * n_mu - 1
        if within_range:
            max_err_within = max(max_err_within, err)
        rows.append((k, numeric, exact, err, within_range))
    return rows, max_err_within


def spherical_harmonic_orthonormality_errors(n_mu, n_phi, azimuthal, n_max):
    quad = sphere_quadrature(n_mu, n_phi, azimuthal=azimuthal)

    worst_err = 0.0
    n_pairs_checked = 0
    err_records = []
    for n1 in range(0, n_max + 1):
        for m1 in range(-n1, n1 + 1):
            for n2 in range(0, n_max + 1):
                for m2 in range(-n2, n2 + 1):
                    def integrand(THETA, PHI, n1=n1, m1=m1, n2=n2, m2=m2):
                        Y1 = sph_harm_y(n1, m1, THETA, PHI)
                        Y2 = sph_harm_y(n2, m2, THETA, PHI)
                        return Y1 * np.conj(Y2)

                    val = integrate_sphere(integrand, quad)
                    is_diagonal = (n1 == n2 and m1 == m2)
                    expected = 1.0 if is_diagonal else 0.0
                    err = abs(val - expected)
                    worst_err = max(worst_err, err)
                    n_pairs_checked += 1
                    err_records.append((n1, m1, n2, m2, err, is_diagonal))

    return worst_err, n_pairs_checked, err_records


def orthonormality_beyond_resolution_error(n_mu, n_phi, azimuthal):
    quad = sphere_quadrature(n_mu, n_phi, azimuthal=azimuthal)
    n_too_high = n_mu + 5  # n1+n2 = 2*n_too_high, well beyond 2*n_mu-1

    curve = []
    for n in range(0, n_too_high + 1):
        def integrand(THETA, PHI, n=n):
            Y = sph_harm_y(n, n, THETA, PHI)
            return Y * np.conj(Y)
        val = integrate_sphere(integrand, quad)
        curve.append((n, abs(val - 1.0)))

    err = curve[-1][1]
    return n_too_high, err, curve


def quadrature_node_view(n_mu, n_phi, azimuthal, color_by,
                          test_function=None, k=None, n=None, m=None):
    quad = sphere_quadrature(n_mu, n_phi, azimuthal=azimuthal)
    THETA, PHI, W = quad["THETA"], quad["PHI"], quad["W"]

    x = (np.sin(THETA) * np.cos(PHI)).ravel()
    y = (np.sin(THETA) * np.sin(PHI)).ravel()
    z = np.cos(THETA).ravel()

    if color_by == "weight":
        color_values = W.ravel()
        color_label = "node weight  w_mu * w_phi"

    elif color_by == "function":
        if test_function == "constant":
            f, f_label = np.ones_like(THETA), "f = 1"
        elif test_function == "mu_power":
            f, f_label = np.cos(THETA) ** k, f"f = mu^{k}"
        elif test_function == "spherical_harmonic":
            f = sph_harm_y(n, m, THETA, PHI).real
            f_label = f"f = Re(Y_{n}^{m})"
        else:
            raise ValueError(f"unknown test_function: '{test_function}'")
        color_values = (W * f).ravel()
        color_label = f"node contribution:  w * ({f_label})"

    else:
        raise ValueError(f"color_by must be 'weight' or 'function', got '{color_by}'")

    return quad, x, y, z, color_values, color_label


def gauss_chebyshev_trapezoidal_equivalence(n_phi):
    phi_trap, w_trap = trapezoidal_azimuthal(n_phi)
    phi_cheb, w_cheb = gauss_chebyshev_azimuthal(n_phi)

    node_shift = np.unique(np.round(phi_cheb - phi_trap, 12))
    weights_identical = bool(np.allclose(w_trap, w_cheb))

    k_rows = []
    for k in range(0, n_phi + 2):
        exact = 2 * np.pi if k == 0 else 0.0
        val_trap = np.sum(w_trap * np.exp(1j * k * phi_trap))
        val_cheb = np.sum(w_cheb * np.exp(1j * k * phi_cheb))
        err_trap = abs(val_trap - exact)
        err_cheb = abs(val_cheb - exact)
        within_range = k <= n_phi - 1
        k_rows.append((k, val_trap.real, val_cheb.real, err_trap, err_cheb, within_range))

    return {
        "phi_trap": phi_trap, "w_trap": w_trap,
        "phi_cheb": phi_cheb, "w_cheb": w_cheb,
        "node_shift": node_shift, "weights_identical": weights_identical,
        "k_rows": k_rows,
    }
