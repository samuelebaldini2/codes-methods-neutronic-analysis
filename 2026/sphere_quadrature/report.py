

import numpy as np


def report_polar_polynomial_exactness(n_mu, rows, max_err_within):
    print(f"\n=== Polynomial exactness in mu, Gauss-Legendre n_mu={n_mu} "
          f"(expected exact up to degree {2 * n_mu - 1}) ===")
    for k, numeric, exact, err, within_range in rows:
        flag = "OK (within range)" if within_range else "out of range"
        if not within_range and err > 1e-8:
            flag += "  <-- indeed NOT exact anymore, as expected"
        print(f"  degree k={k:2d}: numeric={numeric: .12f}  "
              f"exact={exact: .12f}  error={err:.3e}  [{flag}]")
    print(f"  --> max error within the guaranteed range: {max_err_within:.3e}")


def report_spherical_harmonic_orthonormality(n_mu, n_phi, azimuthal, n_max,
                                              worst_err, n_pairs_checked):
    print(f"\n=== Spherical harmonic orthonormality, "
          f"n_mu={n_mu}, n_phi={n_phi}, azimuthal rule='{azimuthal}', "
          f"tested up to n_max={n_max} ===")
    print(f"  (n1,m1,n2,m2) pairs checked: {n_pairs_checked}")
    print(f"  max error over all pairs: {worst_err:.3e}")


def report_orthonormality_beyond_resolution(n_mu, n_phi, n_too_high, err):
    print(f"\n=== Demonstration of the resolution limit "
          f"(n_mu={n_mu}, n_phi={n_phi}) ===")
    print(f"  n={n_too_high} (well beyond the resolution guaranteed by "
          f"n_mu={n_mu}): error = {err:.3e}")
    print("  --> as expected, the error is NO LONGER at machine precision here:"
          " the quadrature no longer resolves this degree exactly.")


def report_quadrature_nodes_3d(n_mu, n_phi, azimuthal, quad, color_by, color_label):
    n_nodes = quad["mu"].size * quad["phi"].size
    total_weight = np.sum(quad["W"])
    print(f"\n=== Quadrature nodes on the unit sphere (3D), "
          f"n_mu={n_mu}, n_phi={n_phi}, azimuthal rule='{azimuthal}' ===")
    print(f"  total nodes: {n_nodes}  ({n_mu} polar x {n_phi} azimuthal)")
    print(f"  sum of all weights (should equal the unit sphere's surface "
          f"area 4*pi = {4 * np.pi:.6f}): {total_weight:.6f}  "
          f"(error: {abs(total_weight - 4 * np.pi):.3e})")
    print(f"  coloring nodes by: {color_by}  ({color_label})")


def report_gauss_chebyshev_trapezoidal_equivalence(n_phi, data):
    print(f"\n=== Gauss-Chebyshev <-> trapezoidal equivalence, "
          f"n_phi={n_phi} (bonus) ===")
    print(f"  trapezoidal nodes:      {data['phi_trap']}")
    print(f"  Gauss-Chebyshev nodes:  {data['phi_cheb']}")
    print(f"  node difference (should be constant = pi/n_phi = "
          f"{np.pi / n_phi:.6f}): {data['node_shift']}")
    print(f"  trapezoidal weights:    {data['w_trap']}")
    print(f"  Gauss-Chebyshev weights:{data['w_cheb']}")
    print(f"  weights identical: {data['weights_identical']}")

    print("  Exactness check on e^{i k phi} for both rules "
          "(the magnitude of the error is expected to be identical; the "
          "sign of the aliased value in Gauss-Chebyshev depends on the "
          "parity of k/n_phi):")
    for k, val_trap, val_cheb, err_trap, err_cheb, within_range in data["k_rows"]:
        print(f"    k={k:2d}  val_trap={val_trap:+.3f}  "
              f"val_cheb={val_cheb:+.3f}  "
              f"err_trap={err_trap:.3e}  err_cheb={err_cheb:.3e}"
              f"  [{'exact as expected' if within_range else 'aliasing as expected'}]")