import numpy as np
import pytest

from sph_harm_compat import sph_harm_y
from sphere_quadrature import (
    gauss_legendre_polar, trapezoidal_azimuthal, gauss_chebyshev_azimuthal,
    sphere_quadrature, integrate_sphere,
)


def analytic_mu_power_integral(k):
    return 0.0 if k % 2 == 1 else 2.0 / (k + 1)


@pytest.mark.parametrize("n_mu", [3, 6, 10])
def test_gauss_legendre_exact_within_degree_2n_minus_1(n_mu):
    mu, w_mu, _ = gauss_legendre_polar(n_mu)
    for k in range(0, 2 * n_mu):  
        numeric = np.sum(w_mu * mu ** k)
        exact = analytic_mu_power_integral(k)
        assert numeric == pytest.approx(exact, abs=1e-12)


@pytest.mark.parametrize("n_mu", [3, 6])
def test_gauss_legendre_not_exact_beyond_degree_2n_minus_1(n_mu):
    mu, w_mu, _ = gauss_legendre_polar(n_mu)
    k = 2 * n_mu  # one degree beyond the guaranteed range
    numeric = np.sum(w_mu * mu ** k)
    exact = analytic_mu_power_integral(k)
    assert abs(numeric - exact) > 1e-6


@pytest.mark.parametrize("azimuthal", ["trapezoidal", "gauss_chebyshev"])
@pytest.mark.parametrize("n_phi", [5, 8, 13])
def test_azimuthal_rules_exact_within_degree_nphi_minus_1(azimuthal, n_phi):
    rule = trapezoidal_azimuthal if azimuthal == "trapezoidal" else gauss_chebyshev_azimuthal
    phi, w_phi = rule(n_phi)
    for k in range(0, n_phi):  
        exact = 2 * np.pi if k == 0 else 0.0
        numeric = np.sum(w_phi * np.exp(1j * k * phi))
        assert numeric == pytest.approx(exact, abs=1e-10)


@pytest.mark.parametrize("n_phi,j", [(8, 1), (8, 2), (8, 3)])
def test_trapezoidal_always_aliases_to_plus_2pi(n_phi, j):
    phi, w_phi = trapezoidal_azimuthal(n_phi)
    k = j * n_phi
    numeric = np.sum(w_phi * np.exp(1j * k * phi))
    assert numeric == pytest.approx(2 * np.pi, abs=1e-10)


@pytest.mark.parametrize("n_phi,j", [(8, 1), (8, 2), (8, 3)])
def test_gauss_chebyshev_aliases_with_sign_minus_one_to_the_j(n_phi, j):
    phi, w_phi = gauss_chebyshev_azimuthal(n_phi)
    k = j * n_phi
    numeric = np.sum(w_phi * np.exp(1j * k * phi))
    expected = ((-1) ** j) * 2 * np.pi
    assert numeric == pytest.approx(expected, abs=1e-10)


def test_gauss_chebyshev_equals_trapezoidal_up_to_phase_shift():
    n_phi = 11
    phi_trap, w_trap = trapezoidal_azimuthal(n_phi)
    phi_cheb, w_cheb = gauss_chebyshev_azimuthal(n_phi)
    assert w_trap == pytest.approx(w_cheb)
    shift = phi_cheb - phi_trap
    assert np.allclose(shift, np.pi / n_phi)


@pytest.mark.parametrize("azimuthal", ["trapezoidal", "gauss_chebyshev"])
def test_spherical_harmonics_orthonormal_within_resolution(azimuthal):
    n_mu, n_phi, n_max = 6, 13, 4
    quad = sphere_quadrature(n_mu, n_phi, azimuthal=azimuthal)

    for n1 in range(n_max + 1):
        for m1 in range(-n1, n1 + 1):
            for n2 in range(n_max + 1):
                for m2 in range(-n2, n2 + 1):
                    def integrand(THETA, PHI, n1=n1, m1=m1, n2=n2, m2=m2):
                        return sph_harm_y(n1, m1, THETA, PHI) * np.conj(
                            sph_harm_y(n2, m2, THETA, PHI)
                        )
                    val = integrate_sphere(integrand, quad)
                    expected = 1.0 if (n1 == n2 and m1 == m2) else 0.0
                    assert val == pytest.approx(expected, abs=1e-10)


def test_spherical_harmonics_not_exact_beyond_resolution():
    n_mu = 6
    quad = sphere_quadrature(n_mu, n_phi=25, azimuthal="trapezoidal")
    n_too_high = n_mu + 5  

    def integrand(THETA, PHI):
        Y = sph_harm_y(n_too_high, n_too_high, THETA, PHI)
        return Y * np.conj(Y)

    val = integrate_sphere(integrand, quad)
    assert abs(val - 1.0) > 1e-3
