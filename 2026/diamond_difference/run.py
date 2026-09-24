
import os
import sys
import numpy as np

import studies
import plotting
from input_reader import read_config, ConfigError


"""Writes everything to several streams at once (console + log file). Like the bash command"""
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

print("Problem definition (from input.txt):")
print(f"  thickness = {cfg.d} cm, cells = {cfg.I}")
print(f"  material  = {cfg.material_mode} "
      f"(sigma_t={cfg.sigma_t}, "
      + (f"c={cfg.c})" if cfg.material_mode == "homogeneous"
         else f"c_left={cfg.c_left}, c_right={cfg.c_right})"))
print(f"  source    = {cfg.source_type} (intensity={cfg.sext}"
      + (f", width={cfg.peak_width} cm)" if cfg.source_type == "peaked" else ")"))
print(f"  BC        = left={cfg.bc_left}, right={cfg.bc_right}")
print(f"  N_quad    = {cfg.N_quad}")
if cfg.g != 0.0 or cfg.source_angular_shape != "isotropic":
    print(f"  anisotropy: g={cfg.g}, source_angular_shape={cfg.source_angular_shape}")
print()

#
#
# 
if cfg.do_solve_and_plot or cfg.do_compare_analytical:
    print("=" * 70)
    print("SOLVE AND PLOT" + (" (comparing boundary conditions)" if cfg.compare_bc else ""))
    print("=" * 70)

    if cfg.compare_bc:
        curves = []
        for bcl, bcr in cfg.bc_compare_list:
            mesh, solver, phi, nit = studies.solve_problem(cfg, cfg.I, bc_left=bcl, bc_right=bcr)
            print(f"  BC={bcl}/{bcr:8s}  iterations={nit:5d}  "
                  f"phi_max={phi.max():9.4f}  phi_min={phi.min():9.4f}")
            curves.append((f"{bcl}/{bcr} ({nit} it.)", mesh, phi))
        title = f"{cfg.material_mode} medium, {cfg.source_type} source\nI={cfg.I} cells"
        out_path = plotting.plot_bc_comparison(cfg.out_dir, curves, title)

    else:
        mesh, solver, phi, nit = studies.solve_problem(cfg, cfg.I)
        print(f"  cells={cfg.I}  iterations={nit}  phi_max={phi.max():.4f}  phi_min={phi.min():.4f}")

        phi0_analytic = None
        if cfg.do_compare_analytical:
            phi0_analytic = cfg.sext / (cfg.sigma_t * (1 - cfg.c))
            err = np.max(np.abs(phi - phi0_analytic))
            print(f"  analytical phi0 = {phi0_analytic:.8f}   max |phi - phi0| = {err:.3e}")

            if cfg.analytical_cell_counts:
                print()
                print("  Convergence check across several mesh sizes:")
                rows = studies.run_analytical_convergence(cfg, phi0_analytic)
                print()

                table_lines = [
                    "Convergence check vs analytical solution",
                    f"analytical phi0 = {phi0_analytic:.8f}",
                    "",
                    f"{'I':>7}  {'iterations':>10}  {'phi':>14}  {'abs_err':>12}",
                ]
                for I_check, nit_check, phi_center, err_check in rows:
                    table_lines.append(
                        f"{I_check:7d}  {nit_check:10d}  {phi_center:14.8f}  {err_check:12.3e}")
                table_path = f"{cfg.out_dir}/analytical_convergence_table.txt"
                with open(table_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(table_lines) + "\n")
                print(f"  -> saved {table_path}")

        title = f"{cfg.material_mode} medium, {cfg.source_type} source, " \
                f"BC={cfg.bc_left}/{cfg.bc_right}\nI={cfg.I} cells"
        out_path = plotting.plot_solution(cfg.out_dir, mesh, phi, nit, title, phi0_analytic)

    print(f"  -> saved {out_path}")
    print()

# 
# 
# 
if cfg.do_grid_convergence:
    print("=" * 70)
    print("GRID CONVERGENCE STUDY")
    print("=" * 70)

    min_edges, plot_results = studies.run_grid_convergence_study(cfg)
    print()

    out_path = plotting.plot_grid_convergence(cfg.out_dir, cfg.gc_cell_counts, min_edges, plot_results)
    print(f"  -> saved {out_path}")
    print()

# 
# 
# 
if cfg.do_iterations_vs_c:
    print("=" * 70)
    print(f"ITERATIONS vs SCATTERING RATIO c ({cfg.ic_c_min} -> {cfg.ic_c_max})"
          + (" (comparing boundary conditions)" if cfg.compare_bc else ""))
    print("=" * 70)

    bc_list = cfg.bc_compare_list if cfg.compare_bc else [(cfg.bc_left, cfg.bc_right)]

    curves = []
    for bcl, bcr in bc_list:
        c_values, iters = studies.run_iterations_vs_c_study(cfg, bcl, bcr)
        print(f"  BC={bcl}/{bcr:8s}  iters(c={c_values[0]})={iters[0]:5d}   "
              f"iters(c={c_values[-1]})={iters[-1]:5d}")
        label = f"{bcl}/{bcr}" if cfg.compare_bc else None
        curves.append((label, c_values, iters))

    title = f"Source-iterations vs scattering ratio\n(I={cfg.I} cells)" if cfg.compare_bc else \
        f"Source-iterations vs scattering ratio\n(BC={cfg.bc_left}/{cfg.bc_right}, I={cfg.I} cells)"
    out_path = plotting.plot_iterations_vs_c(cfg.out_dir, curves, cfg.tol, title, cfg.compare_bc)
    print(f"  -> saved {out_path}")
    print()

#
# 
# 
if cfg.do_plot_anisotropic:
    print("=" * 70)
    print("ANISOTROPIC SCATTERING / SOURCE SHAPE FUNCTIONS")
    print("=" * 70)

    out_path = plotting.plot_anisotropic_functions(
        cfg.out_dir, cfg.g, cfg.sext, cfg.source_angular_shape)
    print(f"  g = {cfg.g}, source_angular_shape = {cfg.source_angular_shape}")
    print(f"  -> saved {out_path}")
    print()

print("Done. Output written to:", os.path.abspath(cfg.out_dir))
print(f"(full console log saved to {_log_path})")
sys.stdout = sys.__stdout__  
_log_file.close()
