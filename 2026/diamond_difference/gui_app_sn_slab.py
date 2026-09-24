"""
Interactive Tkinter GUI for the 1D slab SN solver normally driven by
hand-editing input.txt and running run.py.

Run it with:
    python gui_app.py
"""

from __future__ import annotations

import io
import os
import sys
import types
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import studies
import plotting
from anisotropic import phase_function, source_angular_weight
from input_reader import ConfigError, _parse_int_list, _parse_bc_pairs

_PALETTE = ["tab:blue", "tab:red", "tab:green", "tab:orange", "tab:purple", "tab:brown"]



def draw_anisotropic_functions(fig, g, sext, source_shape):
    mu = np.linspace(-1.0, 1.0, 400)

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(mu, phase_function(mu, g), lw=2, color="tab:blue", label=f"g={g}")
    ax1.axhline(0.5, ls="--", color="gray", lw=1, label="isotropic (g=0)")
    ax1.set_xlabel(r"$\mu_0$")
    ax1.set_ylabel(r"$p(\mu_0)$")
    ax1.set_title("Scattering phase function")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(1, 2, 2)
    s_ext = sext * np.array([source_angular_weight(m, source_shape) for m in mu])
    ax2.plot(mu, s_ext, lw=2, color="tab:red", label=source_shape)
    ax2.axhline(0.5 * sext, ls="--", color="gray", lw=1, label="isotropic")
    ax2.set_xlabel(r"$\mu$")
    ax2.set_ylabel(r"$s_{ext}(\mu)$")
    ax2.set_title("External source angular shape")
    ax2.legend()
    ax2.grid(alpha=0.3)


def draw_bc_comparison(fig, curves, title):
    ax = fig.add_subplot(111)
    for i, (label, mesh, phi) in enumerate(curves):
        ax.plot(mesh.centers, phi, lw=2, color=_PALETTE[i % len(_PALETTE)], label=label)
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(r"$\phi(x)$")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)


def draw_solution(fig, mesh, phi, nit, title, phi0_analytic=None):
    ax = fig.add_subplot(111)
    ax.plot(mesh.centers, phi, lw=2, color="tab:blue", label=f"SN solution ({nit} it.)")
    if phi0_analytic is not None:
        ax.axhline(phi0_analytic, ls="--", color="k",
                   label=f"analytical $\\phi_0={phi0_analytic:.3f}$")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(r"$\phi(x)$")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)


def draw_grid_convergence(fig, gc_cell_counts, min_edges, plot_results):
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.semilogx(gc_cell_counts, min_edges, "o-", color="tab:purple")
    ax1.axhline(0, color="k", lw=1)
    ax1.set_xlabel("number of cells I")
    ax1.set_ylabel("min. edge angular flux over mesh & directions")
    ax1.set_title("Diamond-difference positivity vs mesh refinement")
    ax1.grid(alpha=0.3)

    ax2 = fig.add_subplot(1, 2, 2)
    for I, (mesh, phi) in plot_results.items():
        ax2.plot(mesh.centers, phi, label=f"I={I}")
    ax2.set_xlabel("x [cm]")
    ax2.set_ylabel(r"$\phi(x)$")
    ax2.set_title("Scalar flux convergence with mesh refinement")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)


def draw_iterations_vs_c(fig, curves, tol, title, show_legend):
    ax = fig.add_subplot(111)
    for i, (label, c_values, iters) in enumerate(curves):
        ax.semilogy(c_values, iters, "o-", ms=3, color=_PALETTE[i % len(_PALETTE)], label=label)
    ax.set_xlabel("scattering ratio c")
    ax.set_ylabel(f"number of source-iterations to reach tol={tol:.0e}")
    ax.set_title(title)
    if show_legend:
        ax.legend()
    ax.grid(alpha=0.3, which="both")


class SlabSNApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("1D Slab SN Solver \u2014 Interactive GUI")
        self.root.geometry("1340x860")
        self.root.minsize(1020, 660)

        self.field_vars: dict[str, tk.Variable] = {}
        self._last_plots = None       # list of {"title","draw","save"} from the last run
        self._last_analytical_table = None
        self._param_panel_visible = True

        self._build_layout()
        self._build_param_panel()
        self._update_material_visibility()
        self._update_source_visibility()
        self._update_action_visibility()
        self._run()  # populate with sensible defaults on startup


    def _build_layout(self):
        bar = ttk.Frame(self.root, padding=8)
        bar.pack(side=tk.TOP, fill=tk.X)

        self.run_button = ttk.Button(bar, text="Run", command=self._run)
        self.run_button.pack(side=tk.LEFT)

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)

        ttk.Label(bar, text="Output folder:").pack(side=tk.LEFT)
        self.out_dir_var = tk.StringVar(value="out")
        ttk.Entry(bar, textvariable=self.out_dir_var, width=22).pack(side=tk.LEFT, padx=(4, 4))
        ttk.Button(bar, text="Browse...", command=self._browse_out_dir).pack(side=tk.LEFT)
        ttk.Button(bar, text="Save results...", command=self._save_results).pack(side=tk.LEFT, padx=(12, 0))

        self.toggle_param_button = ttk.Button(
            bar, text="\u25c0 Hide parameters", command=self._toggle_param_panel)
        self.toggle_param_button.pack(side=tk.LEFT, padx=(12, 0))

        self.body = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.body.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        # --- left: scrollable parameter panel -----------------------------
        self.left_frame = ttk.Frame(self.body)
        self.body.add(self.left_frame, weight=0)
        left = self.left_frame

        self.param_canvas = tk.Canvas(left, borderwidth=0, highlightthickness=0, width=340)
        vsb = ttk.Scrollbar(left, orient="vertical", command=self.param_canvas.yview)
        self.param_inner = ttk.Frame(self.param_canvas)
        self.param_inner.bind(
            "<Configure>",
            lambda e: self.param_canvas.configure(scrollregion=self.param_canvas.bbox("all")))
        self.param_canvas.create_window((0, 0), window=self.param_inner, anchor="nw")
        self.param_canvas.configure(yscrollcommand=vsb.set)
        self.param_canvas.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.param_inner.columnconfigure(0, weight=1)

        def _wheel(event):
            delta = -1 * (event.delta // 120) if event.delta else (-1 if event.num == 4 else 1)
            self.param_canvas.yview_scroll(int(delta), "units")
        self.param_canvas.bind("<MouseWheel>", _wheel)
        self.param_canvas.bind("<Button-4>", _wheel)
        self.param_canvas.bind("<Button-5>", _wheel)

        # --- right: notebook (plot tabs + Report)
        right = ttk.Frame(self.body)
        self.body.add(right, weight=1)

        self.notebook = ttk.Notebook(right)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.report_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.report_frame, text="Report")
        self.report_text = tk.Text(self.report_frame, wrap="word", state="disabled", font=("Consolas", 10))
        rvsb = ttk.Scrollbar(self.report_frame, orient="vertical", command=self.report_text.yview)
        self.report_text.configure(yscrollcommand=rvsb.set)
        self.report_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4, 0), pady=4)
        rvsb.pack(side=tk.RIGHT, fill=tk.Y, pady=4)

        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W,
                  padding=(6, 2)).pack(side=tk.BOTTOM, fill=tk.X)

        # give the parameter pane a sensible starting width (drag to resize)
        self.root.after(0, lambda: self.body.sashpos(0, 380))

    def _toggle_param_panel(self):
        if self._param_panel_visible:
            self.body.forget(self.left_frame)
            self.toggle_param_button.config(text="\u25b6 Show parameters")
        else:
            self.body.insert(0, self.left_frame, weight=0)
            self.root.after(0, lambda: self.body.sashpos(0, 380))
            self.toggle_param_button.config(text="\u25c0 Hide parameters")
        self._param_panel_visible = not self._param_panel_visible

    def _section(self, row, title):
        frame = ttk.LabelFrame(self.param_inner, text=title, padding=8)
        frame.grid(row=row, column=0, sticky="ew", pady=(0, 8), padx=2)
        frame.columnconfigure(1, weight=1)
        return frame

    def _label(self, parent, row, text, columnspan=1):
        ttk.Label(parent, text=text, wraplength=210, justify="left").grid(
            row=row, column=0, columnspan=columnspan, sticky="w", pady=(4, 0))

    def _entry(self, parent, row, key, label, default, width=14):
        self._label(parent, row, label)
        var = tk.StringVar(value=str(default))
        ttk.Entry(parent, textvariable=var, width=width).grid(
            row=row, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))
        self.field_vars[key] = var
        return row + 1

    def _choice(self, parent, row, key, label, choices, default, on_change=None, width=16):
        self._label(parent, row, label)
        var = tk.StringVar(value=default)
        combo = ttk.Combobox(parent, textvariable=var, values=list(choices),
                              state="readonly", width=width)
        combo.grid(row=row, column=1, sticky="ew", padx=(6, 0), pady=(4, 0))
        if on_change is not None:
            combo.bind("<<ComboboxSelected>>", lambda e: on_change())
        self.field_vars[key] = var
        return row + 1

    def _check(self, parent, row, key, label, default, on_change=None):
        var = tk.BooleanVar(value=default)
        cb = ttk.Checkbutton(parent, text=label, variable=var,
                              command=on_change if on_change else None)
        cb.grid(row=row, column=0, columnspan=2, sticky="w", pady=(4, 0))
        self.field_vars[key] = var
        return row + 1

    @staticmethod
    def _set_visible(frame, visible):
        if visible:
            frame.grid()
        else:
            frame.grid_remove()

    def _build_param_panel(self):
        srow = 0

        # GEOMETRY 
        f = self._section(srow, "Geometry"); srow += 1
        row = 0
        row = self._entry(f, row, "thickness_cm", "thickness_cm  [cm]", 5.0)
        row = self._entry(f, row, "number_of_cells", "number_of_cells", 100)

        # MATERIAL 
        f = self._section(srow, "Material"); srow += 1
        row = 0
        row = self._choice(f, row, "mode", "mode", ("homogeneous", "heterogeneous"),
                            "homogeneous", on_change=self._update_material_visibility)
        row = self._entry(f, row, "sigma_t_cm-1", "sigma_t_cm-1  [cm^-1]", 1.8)

        self.material_homog_frame = ttk.Frame(f)
        self.material_homog_frame.grid(row=row, column=0, columnspan=2, sticky="ew")
        self.material_homog_frame.columnconfigure(1, weight=1)
        self._entry(self.material_homog_frame, 0, "c", "c  (Sigma_s/Sigma_t)", 0.9)
        row += 1

        self.material_heterog_frame = ttk.Frame(f)
        self.material_heterog_frame.grid(row=row, column=0, columnspan=2, sticky="ew")
        self.material_heterog_frame.columnconfigure(1, weight=1)
        r2 = self._entry(self.material_heterog_frame, 0, "c_left", "c_left (x<d/2)", 0.5)
        self._entry(self.material_heterog_frame, r2, "c_right", "c_right (x>d/2)", 0.99)
        row += 1

        row = self._entry(f, row, "g", "g  (mean cosine, -1/3..1/3)", 0.2)

        # SOURCE 
        f = self._section(srow, "Source"); srow += 1
        row = 0
        row = self._choice(f, row, "type", "type", ("uniform", "peaked"),
                            "uniform", on_change=self._update_source_visibility)
        row = self._entry(f, row, "intensity", "intensity  Q [cm^-3 s^-1]", 10.0)

        self.source_peaked_frame = ttk.Frame(f)
        self.source_peaked_frame.grid(row=row, column=0, columnspan=2, sticky="ew")
        self.source_peaked_frame.columnconfigure(1, weight=1)
        self._entry(self.source_peaked_frame, 0, "peaked_width_cm", "peaked_width_cm", 0.1)
        row += 1

        row = self._choice(f, row, "angular_shape", "angular_shape",
                            ("isotropic", "forward_peaked"), "forward_peaked")

        # BOUNDARY CONDITIONS 
        f = self._section(srow, "Boundary conditions"); srow += 1
        row = 0
        row = self._choice(f, row, "left", "left", ("void", "reflect"), "void")
        row = self._choice(f, row, "right", "right", ("void", "reflect"), "void")
        row = self._entry(f, row, "compare_list", "compare_list (only used if 'compare "
                           "boundary conditions' is checked below; e.g. "
                           "'void-void, reflect-reflect')",
                           "void-void, reflect-reflect, void-reflect", width=26)

        # ANGULAR QUADRATURE 
        f = self._section(srow, "Angular quadrature"); srow += 1
        self._entry(f, 0, "order_N", "order_N  (even)", 4)

        # SOLVER 
        f = self._section(srow, "Solver"); srow += 1
        row = 0
        row = self._entry(f, row, "tolerance", "tolerance", "1e-8")
        row = self._entry(f, row, "max_iterations", "max_iterations", 20000)

        # ACTIONS 
        f = self._section(srow, "Actions"); srow += 1
        row = 0
        row = self._check(f, row, "solve_and_plot", "Solve and plot", True)
        row = self._check(f, row, "compare_with_analytical", "Compare with analytical",
                           False, on_change=self._update_action_visibility)
        row = self._check(f, row, "grid_convergence_study", "Grid convergence study",
                           False, on_change=self._update_action_visibility)
        row = self._check(f, row, "iterations_vs_c_study", "Iterations vs c study",
                           False, on_change=self._update_action_visibility)
        row = self._check(f, row, "compare_boundary_conditions",
                           "Compare boundary conditions (overlay compare_list above)", True)
        row = self._check(f, row, "plot_anisotropic_functions",
                           "Plot anisotropic functions", True)

        # ANALYTICAL_COMPARISON  
        self.analytical_frame = self._section(srow, "Analytical comparison"); srow += 1
        self._entry(self.analytical_frame, 0, "extra_cell_counts",
                     "extra_cell_counts (comma-separated, optional)", "50, 100, 200, 400", width=24)

        # GRID_CONVERGENCE_STUDY 
        self.grid_conv_frame = self._section(srow, "Grid convergence study"); srow += 1
        row = 0
        row = self._entry(self.grid_conv_frame, row, "cell_counts",
                           "cell_counts (comma-separated)",
                           "3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 80, 120, 200, 300", width=24)
        self._entry(self.grid_conv_frame, row, "cell_counts_for_plot",
                    "cell_counts_for_plot", "4, 10, 20, 50, 200", width=24)

        # ITERATIONS_VS_C_STUDY
        self.iter_c_frame = self._section(srow, "Iterations vs c study"); srow += 1
        row = 0
        row = self._entry(self.iter_c_frame, row, "c_min", "c_min", 0.20)
        row = self._entry(self.iter_c_frame, row, "c_threshold", "c_threshold", 0.90)
        row = self._entry(self.iter_c_frame, row, "c_max", "c_max", 0.99)
        row = self._entry(self.iter_c_frame, row, "coarse_step", "coarse_step", 0.05)
        row = self._entry(self.iter_c_frame, row, "fine_step", "fine_step", 0.01)
        self._entry(self.iter_c_frame, row, "ic_max_iterations", "max_iterations", 50000)

    def _update_material_visibility(self):
        homog = self.field_vars["mode"].get() == "homogeneous"
        self._set_visible(self.material_homog_frame, homog)
        self._set_visible(self.material_heterog_frame, not homog)

    def _update_source_visibility(self):
        self._set_visible(self.source_peaked_frame, self.field_vars["type"].get() == "peaked")

    def _update_action_visibility(self):
        self._set_visible(self.analytical_frame, self.field_vars["compare_with_analytical"].get())
        self._set_visible(self.grid_conv_frame, self.field_vars["grid_convergence_study"].get())
        self._set_visible(self.iter_c_frame, self.field_vars["iterations_vs_c_study"].get())


    def _float(self, key, label=None):
        raw = self.field_vars[key].get().strip()
        try:
            return float(raw)
        except ValueError:
            raise ValueError(f"'{label or key}' must be a number, found: '{raw}'")

    def _int(self, key, label=None):
        raw = self.field_vars[key].get().strip()
        try:
            return int(raw)
        except ValueError:
            raise ValueError(f"'{label or key}' must be an integer, found: '{raw}'")

    def _collect_config(self):
        cfg = types.SimpleNamespace()

        # GEOMETRY
        cfg.d = self._float("thickness_cm")
        if cfg.d <= 0:
            raise ValueError("'thickness_cm' must be positive")
        cfg.I = self._int("number_of_cells")
        if cfg.I <= 0:
            raise ValueError("'number_of_cells' must be a positive integer")

        # MATERIAL
        cfg.material_mode = self.field_vars["mode"].get()
        cfg.sigma_t = self._float("sigma_t_cm-1")
        if cfg.sigma_t <= 0:
            raise ValueError("'sigma_t_cm-1' must be positive")
        if cfg.material_mode == "homogeneous":
            cfg.c = self._float("c")
            if not (0.0 <= cfg.c < 1.0):
                raise ValueError("'c' must be between 0 (included) and 1 (excluded)")
            cfg.c_left = cfg.c_right = cfg.c
        else:
            cfg.c_left = self._float("c_left")
            cfg.c_right = self._float("c_right")
            for name, val in [("c_left", cfg.c_left), ("c_right", cfg.c_right)]:
                if not (0.0 <= val < 1.0):
                    raise ValueError(f"'{name}' must be between 0 (included) and 1 (excluded)")
            cfg.c = None
        cfg.g = self._float("g")
        if not (-1.0 / 3.0 <= cfg.g <= 1.0 / 3.0):
            raise ValueError(
                "'g' must satisfy -1/3 <= g <= 1/3 (otherwise the scattering phase "
                "function p(mu0)=1/2(1+3 g mu0) would be negative for some mu0)")

        # SOURCE
        cfg.source_type = self.field_vars["type"].get()
        cfg.sext = self._float("intensity")
        cfg.source_angular_shape = self.field_vars["angular_shape"].get()
        if cfg.source_type == "peaked":
            cfg.peak_width = self._float("peaked_width_cm")
            if not (0 < cfg.peak_width <= cfg.d):
                raise ValueError("'peaked_width_cm' must be between 0 and the slab thickness")
        else:
            cfg.peak_width = None

        # BOUNDARY CONDITIONS
        cfg.bc_left = self.field_vars["left"].get()
        cfg.bc_right = self.field_vars["right"].get()

        # ANGULAR QUADRATURE
        cfg.N_quad = self._int("order_N")
        if cfg.N_quad <= 0 or cfg.N_quad % 2 != 0:
            raise ValueError("'order_N' must be a positive even integer")

        # SOLVER
        cfg.tol = self._float("tolerance")
        cfg.max_iter = self._int("max_iterations")

        # ACTIONS
        cfg.do_solve_and_plot = self.field_vars["solve_and_plot"].get()
        cfg.do_compare_analytical = self.field_vars["compare_with_analytical"].get()
        cfg.do_grid_convergence = self.field_vars["grid_convergence_study"].get()
        cfg.do_iterations_vs_c = self.field_vars["iterations_vs_c_study"].get()
        cfg.compare_bc = self.field_vars["compare_boundary_conditions"].get()
        cfg.do_plot_anisotropic = self.field_vars["plot_anisotropic_functions"].get()

        if not any([cfg.do_solve_and_plot, cfg.do_compare_analytical,
                    cfg.do_grid_convergence, cfg.do_iterations_vs_c,
                    cfg.do_plot_anisotropic]):
            raise ValueError("All actions are unchecked: nothing to do. Enable at least one.")

        if cfg.compare_bc:
            if cfg.do_compare_analytical:
                raise ValueError(
                    "'Compare boundary conditions' cannot be combined with "
                    "'Compare with analytical' (the analytical check assumes a single, "
                    "fixed reflect/reflect configuration). Turn one of the two off.")
            if not (cfg.do_solve_and_plot or cfg.do_iterations_vs_c):
                raise ValueError(
                    "'Compare boundary conditions' has no effect unless 'Solve and plot' "
                    "and/or 'Iterations vs c study' is also enabled.")
            raw_bc_list = self.field_vars["compare_list"].get()
            cfg.bc_compare_list = _parse_bc_pairs(raw_bc_list, "BOUNDARY_CONDITIONS", "compare_list")

        if cfg.do_compare_analytical:
            if cfg.material_mode != "homogeneous":
                raise ValueError("'Compare with analytical' requires material mode = homogeneous")
            if cfg.bc_left != "reflect" or cfg.bc_right != "reflect":
                raise ValueError(
                    "'Compare with analytical' requires both boundary conditions to be "
                    "'reflect' (reflect/reflect approximates an infinite homogeneous medium)")
            if cfg.g != 0.0 or cfg.source_angular_shape != "isotropic":
                raise ValueError(
                    "'Compare with analytical' requires isotropic scattering and source "
                    "(set g = 0.0 and angular_shape = isotropic, or turn this action off)")
            raw_extra_I = self.field_vars["extra_cell_counts"].get().strip()
            cfg.analytical_cell_counts = (
                _parse_int_list(raw_extra_I, "ANALYTICAL_COMPARISON", "extra_cell_counts")
                if raw_extra_I else None)

        if cfg.do_grid_convergence:
            raw_Ilist = self.field_vars["cell_counts"].get()
            cfg.gc_cell_counts = _parse_int_list(raw_Ilist, "GRID_CONVERGENCE_STUDY", "cell_counts")
            raw_Iplot = self.field_vars["cell_counts_for_plot"].get()
            cfg.gc_cell_counts_plot = _parse_int_list(
                raw_Iplot, "GRID_CONVERGENCE_STUDY", "cell_counts_for_plot")

        if cfg.do_iterations_vs_c:
            cfg.ic_c_min = self._float("c_min")
            cfg.ic_c_threshold = self._float("c_threshold")
            cfg.ic_c_max = self._float("c_max")
            cfg.ic_step_coarse = self._float("coarse_step")
            cfg.ic_step_fine = self._float("fine_step")
            cfg.ic_max_iter = self._int("ic_max_iterations", label="max_iterations")

        cfg.out_dir = self.out_dir_var.get().strip() or "out"
        return cfg

    def _run(self):
        try:
            cfg = self._collect_config()
        except (ValueError, ConfigError) as e:
            messagebox.showerror("Invalid input", str(e))
            return

        self.run_button.config(state=tk.DISABLED)
        self.status_var.set("Running...")
        self.root.update_idletasks()

        buf = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buf
        plots_made = []
        self._last_analytical_table = None
        try:
            self._run_pipeline(cfg, plots_made)
            print("(plots are ready in the tabs on the right; "
                  "use 'Save results...' to write them to disk)")
        except Exception as e:  
            sys.stdout = old_stdout
            self.run_button.config(state=tk.NORMAL)
            messagebox.showerror("Run failed", str(e))
            self.status_var.set("Error.")
            return
        finally:
            sys.stdout = old_stdout

        self._rebuild_plot_tabs(plots_made)
        self._set_report_text(buf.getvalue())
        self._last_plots = plots_made
        self.run_button.config(state=tk.NORMAL)
        self.status_var.set(f"Done: {len(plots_made)} plot(s) produced.")

    def _run_pipeline(self, cfg, plots_made):
        print("Problem definition:")
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
                plots_made.append({
                    "title": "Solution (BC comparison)",
                    "draw": lambda fig, curves=curves, title=title: draw_bc_comparison(fig, curves, title),
                    "save": lambda out_dir, curves=curves, title=title:
                        plotting.plot_bc_comparison(out_dir, curves, title),
                })
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
                        self._last_analytical_table = (rows, phi0_analytic)

                title = f"{cfg.material_mode} medium, {cfg.source_type} source, " \
                        f"BC={cfg.bc_left}/{cfg.bc_right}\nI={cfg.I} cells"
                plots_made.append({
                    "title": "Solution",
                    "draw": lambda fig, mesh=mesh, phi=phi, nit=nit, title=title, pa=phi0_analytic:
                        draw_solution(fig, mesh, phi, nit, title, pa),
                    "save": lambda out_dir, mesh=mesh, phi=phi, nit=nit, title=title, pa=phi0_analytic:
                        plotting.plot_solution(out_dir, mesh, phi, nit, title, pa),
                })
            print()

        if cfg.do_grid_convergence:
            print("=" * 70)
            print("GRID CONVERGENCE STUDY")
            print("=" * 70)
            min_edges, plot_results = studies.run_grid_convergence_study(cfg)
            print()
            plots_made.append({
                "title": "Grid convergence",
                "draw": lambda fig, a=cfg.gc_cell_counts, b=min_edges, c=plot_results:
                    draw_grid_convergence(fig, a, b, c),
                "save": lambda out_dir, a=cfg.gc_cell_counts, b=min_edges, c=plot_results:
                    plotting.plot_grid_convergence(out_dir, a, b, c),
            })
            print()

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

            title = (f"Source-iterations vs scattering ratio\n(I={cfg.I} cells)" if cfg.compare_bc else
                     f"Source-iterations vs scattering ratio\n(BC={cfg.bc_left}/{cfg.bc_right}, I={cfg.I} cells)")
            plots_made.append({
                "title": "Iterations vs c",
                "draw": lambda fig, c=curves, tol=cfg.tol, t=title, sl=cfg.compare_bc:
                    draw_iterations_vs_c(fig, c, tol, t, sl),
                "save": lambda out_dir, c=curves, tol=cfg.tol, t=title, sl=cfg.compare_bc:
                    plotting.plot_iterations_vs_c(out_dir, c, tol, t, sl),
            })
            print()

        if cfg.do_plot_anisotropic:
            print("=" * 70)
            print("ANISOTROPIC SCATTERING / SOURCE SHAPE FUNCTIONS")
            print("=" * 70)
            print(f"  g = {cfg.g}, source_angular_shape = {cfg.source_angular_shape}")
            plots_made.append({
                "title": "Anisotropic functions",
                "draw": lambda fig, g=cfg.g, s=cfg.sext, sh=cfg.source_angular_shape:
                    draw_anisotropic_functions(fig, g, s, sh),
                "save": lambda out_dir, g=cfg.g, s=cfg.sext, sh=cfg.source_angular_shape:
                    plotting.plot_anisotropic_functions(out_dir, g, s, sh),
            })
            print()


    def _rebuild_plot_tabs(self, plots_made):
        for tab_id in self.notebook.tabs():
            if self.notebook.tab(tab_id, "text") != "Report":
                self.notebook.forget(tab_id)

        for i, item in enumerate(plots_made):
            frame = ttk.Frame(self.notebook)
            self.notebook.insert(i, frame, text=item["title"])
            fig = Figure(figsize=(8, 5.2), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            toolbar_frame = ttk.Frame(frame)
            toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)
            toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
            item["draw"](fig)
            canvas.draw()

        if plots_made:
            self.notebook.select(0)
        else:
            self.notebook.select(self.report_frame)

    def _set_report_text(self, text):
        self.report_text.config(state="normal")
        self.report_text.delete("1.0", tk.END)
        self.report_text.insert(tk.END, text)
        self.report_text.config(state="disabled")

    def _browse_out_dir(self):
        directory = filedialog.askdirectory(title="Choose output folder")
        if directory:
            self.out_dir_var.set(directory)

    def _save_results(self):
        if not self._last_plots:
            messagebox.showinfo("Nothing to save", "Run the solver first.")
            return
        out_dir = self.out_dir_var.get().strip() or "out"
        try:
            os.makedirs(out_dir, exist_ok=True)
            saved_paths = [item["save"](out_dir) for item in self._last_plots]

            log_path = os.path.join(out_dir, "run_log.txt")
            with open(log_path, "w", encoding="utf-8") as f:
                f.write(self.report_text.get("1.0", tk.END))
            saved_paths.append(log_path)

            if self._last_analytical_table:
                rows, phi0_analytic = self._last_analytical_table
                lines = [
                    "Convergence check vs analytical solution",
                    f"analytical phi0 = {phi0_analytic:.8f}",
                    "",
                    f"{'I':>7}  {'iterations':>10}  {'phi':>14}  {'abs_err':>12}",
                ]
                for I_check, nit_check, phi_center, err_check in rows:
                    lines.append(
                        f"{I_check:7d}  {nit_check:10d}  {phi_center:14.8f}  {err_check:12.3e}")
                table_path = os.path.join(out_dir, "analytical_convergence_table.txt")
                with open(table_path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines) + "\n")
                saved_paths.append(table_path)
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Save failed", str(e))
            return

        self.status_var.set(f"Saved {len(saved_paths)} file(s) to {os.path.abspath(out_dir)}")
        messagebox.showinfo("Saved", "Wrote:\n" + "\n".join(saved_paths))


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    SlabSNApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
