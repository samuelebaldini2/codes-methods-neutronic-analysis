import os
import sys
import numpy as np

import checks
import report
import plots
from input_reader import read_config, ConfigError


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)

    def flush(self):
        for s in self.streams:
            s.flush()


try:
    cfg = read_config("input.txt")
except ConfigError as e:
    raise SystemExit(f"ERROR in input file: {e}")

os.makedirs(cfg.out_dir, exist_ok=True)


_log_path = f"{cfg.out_dir}/run_log.txt"
_log_file = open(_log_path, "w", encoding="utf-8")
sys.stdout = _Tee(sys.stdout, _log_file)

np.set_printoptions(precision=3, suppress=False)

print(f"Selected check (from input.txt): {cfg.which}")


_plot_path = f"{cfg.out_dir}/plot_{cfg.which}.png"

if cfg.which == "polar_polynomial_exactness":
    rows, max_err_within = checks.polar_polynomial_errors(n_mu=cfg.n_mu, k_max=cfg.k_max)
    report.report_polar_polynomial_exactness(n_mu=cfg.n_mu, rows=rows, max_err_within=max_err_within)
    plots.plot_polar_polynomial_exactness(n_mu=cfg.n_mu, rows=rows, out_path=_plot_path)
    assert max_err_within < 1e-12, "should have been exact within machine precision!"

elif cfg.which == "spherical_harmonic_orthonormality":
    worst_err, n_pairs, err_records = checks.spherical_harmonic_orthonormality_errors(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal, n_max=cfg.n_max)
    report.report_spherical_harmonic_orthonormality(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal, n_max=cfg.n_max,
        worst_err=worst_err, n_pairs_checked=n_pairs)
    plots.plot_spherical_harmonic_orthonormality(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal, n_max=cfg.n_max,
        err_records=err_records, out_path=_plot_path)
    assert worst_err < 1e-10, "expected an error at machine precision within the resolution"

elif cfg.which == "orthonormality_beyond_resolution":
    n_too_high, err_beyond, curve = checks.orthonormality_beyond_resolution_error(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal)
    report.report_orthonormality_beyond_resolution(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, n_too_high=n_too_high, err=err_beyond)
    plots.plot_orthonormality_beyond_resolution(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, curve=curve, n_too_high=n_too_high,
        err=err_beyond, out_path=_plot_path)
    assert err_beyond > 1e-6, "expected an appreciable error beyond the resolution"

elif cfg.which == "gauss_chebyshev_trapezoidal_equivalence":
    equivalence_data = checks.gauss_chebyshev_trapezoidal_equivalence(n_phi=cfg.n_phi)
    report.report_gauss_chebyshev_trapezoidal_equivalence(n_phi=cfg.n_phi, data=equivalence_data)
    plots.plot_gauss_chebyshev_trapezoidal_equivalence(
        n_phi=cfg.n_phi, data=equivalence_data, out_path=_plot_path)

else:  
    quad, x, y, z, color_values, color_label = checks.quadrature_node_view(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal,
        color_by=cfg.color_by, test_function=cfg.test_function,
        k=cfg.k, n=cfg.n, m=cfg.m)
    report.report_quadrature_nodes_3d(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal,
        quad=quad, color_by=cfg.color_by, color_label=color_label)
    plots.plot_quadrature_nodes_3d(
        n_mu=cfg.n_mu, n_phi=cfg.n_phi, azimuthal=cfg.azimuthal,
        x=x, y=y, z=z, color_values=color_values, color_label=color_label,
        out_path=_plot_path)

print()
print("Check completed successfully.")
print("Done. Output written to:", os.path.abspath(cfg.out_dir))
print(f"(full console log saved to {_log_path})")
print(f"(plot saved to {_plot_path})")
sys.stdout = sys.__stdout__  
_log_file.close()
