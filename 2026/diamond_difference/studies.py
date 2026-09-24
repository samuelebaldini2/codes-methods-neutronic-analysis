#Case studies
import numpy as np

from sn_solver import Mesh, make_uniform_source, make_peaked_source, SNSlabSolver
from anisotropic import AnisotropicSNSlabSolver, direction_source


def _is_anisotropic(cfg):
    return (getattr(cfg, "g", 0.0) != 0.0
            or getattr(cfg, "source_angular_shape", "isotropic") != "isotropic")


def build_mesh_and_source(cfg, I):
    mesh = Mesh(cfg.d, I)
    if cfg.material_mode == "homogeneous":
        c_array = cfg.c
    else:
        c_array = np.where(mesh.centers < cfg.d / 2, cfg.c_left, cfg.c_right)
    if cfg.source_type == "uniform":
        Q = make_uniform_source(mesh, cfg.sext)
    else:
        Q, _actual_w, _ncells = make_peaked_source(mesh, cfg.sext, cfg.peak_width)
    return mesh, c_array, Q

############################
def solve_problem(cfg, I, bc_left=None, bc_right=None, c_override=None, tol=None, max_iter=None):
    mesh, c_array, Q = build_mesh_and_source(cfg, I)
    if c_override is not None:
        c_array = c_override
    bcl = bc_left if bc_left is not None else cfg.bc_left
    bcr = bc_right if bc_right is not None else cfg.bc_right
    if _is_anisotropic(cfg):
        solver = AnisotropicSNSlabSolver(
            mesh, cfg.sigma_t, c_array, Q, N=cfg.N_quad,
            bc_left=bcl, bc_right=bcr,
            g=cfg.g, source_shape=cfg.source_angular_shape,
        )
    else:
        solver = SNSlabSolver(
            mesh, cfg.sigma_t, c_array, Q, N=cfg.N_quad,
            bc_left=bcl, bc_right=bcr,
        )
    phi, nit = solver.solve(
        tol=tol if tol is not None else cfg.tol,
        max_iter=max_iter if max_iter is not None else cfg.max_iter,
    )
    return mesh, solver, phi, nit
############################

def min_edge_flux(solver, phi):
    mesh = solver.mesh
    I = mesh.I
    dx = mesh.dx
    sigma_t = solver.sigma_t
    is_aniso = hasattr(solver, "sigma_s1")  # True only for AnisotropicSNSlabSolver
    phi1 = solver.phi1 if is_aniso else None

    def _Q(mu_signed):
        if is_aniso:
            return direction_source(mu_signed, phi, phi1, solver.sigma_s,
                                     solver.sigma_s1, solver.Q_ext, solver.source_shape)
        return 0.5 * (solver.sigma_s * phi + solver.Q_ext)

    min_edge = np.inf
    for m in solver.pos:
        mu_m = solver.mu[m]
        Q_dir = _Q(mu_m)
        psi_edge = 0.0
        coef = 2.0 * mu_m / dx
        for i in range(I):
            denom = coef[i] + sigma_t[i]
            psi_i = (Q_dir[i] + coef[i] * psi_edge) / denom
            psi_edge = 2 * psi_i - psi_edge
            min_edge = min(min_edge, psi_edge)
    for m in solver.neg:
        mu_m = solver.mu[m]
        mu_abs = -mu_m
        Q_dir = _Q(mu_m)
        psi_edge = 0.0
        coef = 2.0 * mu_abs / dx
        for i in range(I - 1, -1, -1):
            denom = coef[i] + sigma_t[i]
            psi_i = (Q_dir[i] + coef[i] * psi_edge) / denom
            psi_edge = 2 * psi_i - psi_edge
            min_edge = min(min_edge, psi_edge)
    return min_edge


"""Runs the optional multi-mesh convergence check against the
analytical infinite-medium solution, printing one progress line
per mesh size. Returns a list of (I, n_iter, phi_center, abs_err)
tuples (empty list if cfg.analytical_cell_counts is not set)."""
def run_analytical_convergence(cfg, phi0_analytic):
    rows = []
    if not cfg.analytical_cell_counts:
        return rows
    for I_check in cfg.analytical_cell_counts:
        _m, _s, phi_check, nit_check = solve_problem(cfg, I_check)
        err_check = np.max(np.abs(phi_check - phi0_analytic))
        rows.append((I_check, nit_check, phi_check[0], err_check))
        print(f"    I={I_check:5d}  iterations={nit_check:5d}  "
              f"phi={phi_check[0]:.8f}  |err|={err_check:.3e}")
    return rows

"""Solves the problem on every mesh size in cfg.gc_cell_counts and
tracks the minimum edge angular flux on each, to expose the
diamond-difference negativity artifact as the mesh is refined.
Also solves on cfg.gc_cell_counts_plot for the flux-overlay plot."""
def run_grid_convergence_study(cfg):
    min_edges = []
    for I in cfg.gc_cell_counts:
        mesh, solver, phi, nit = solve_problem(cfg, I)
        me = min_edge_flux(solver, phi)
        min_edges.append(me)
        flag = " <-- NEGATIVE EDGE FLUX" if me < 0 else ""
        print(f"  I={I:4d}  dx={cfg.d/I:.4f} cm  min edge psi = {me: .4e}{flag}")

    plot_results = {}
    for I in cfg.gc_cell_counts_plot:
        mesh, solver, phi, nit = solve_problem(cfg, I)
        plot_results[I] = (mesh, phi)

    return min_edges, plot_results


"""Scans the scattering ratio c (coarse step below cfg.ic_c_threshold,
fine step above it) for one boundary-condition pair.
"""
def run_iterations_vs_c_study(cfg, bc_left, bc_right):
    c_values = np.concatenate([
        np.arange(cfg.ic_c_min, cfg.ic_c_threshold, cfg.ic_step_coarse),
        np.arange(cfg.ic_c_threshold, cfg.ic_c_max + cfg.ic_step_fine / 2, cfg.ic_step_fine),
    ])
    c_values = np.unique(np.round(c_values, 4))

    iters = []
    for c in c_values:
        c_array = np.full(cfg.I, c) if cfg.material_mode == "heterogeneous" else c
        _mesh, _solver, _phi, nit = solve_problem(
            cfg, cfg.I, bc_left=bc_left, bc_right=bc_right,
            c_override=c_array, max_iter=cfg.ic_max_iter)
        iters.append(nit)
    return c_values, iters
