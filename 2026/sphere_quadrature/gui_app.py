"""
Interactive Tkinter GUI for the sphere-quadrature checks that are normally
driven by editing input.txt and running run.py.

Run it with:
    python gui_app.py
"""

from __future__ import annotations

import io
import os
import sys
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- registers the 3d projection

import checks
import report
import plots
from input_reader import (
    _VALID_CHECKS,
    _VALID_AZIMUTHAL,
    _VALID_COLOR_BY,
    _VALID_TEST_FUNCTION,
)

CHECK_LABELS = {
    "polar_polynomial_exactness": "Polar polynomial exactness",
    "spherical_harmonic_orthonormality": "Spherical harmonic orthonormality",
    "orthonormality_beyond_resolution": "Orthonormality beyond resolution",
    "gauss_chebyshev_trapezoidal_equivalence": "Gauss-Chebyshev vs trapezoidal",
    "quadrature_nodes_3d": "Quadrature nodes (3D view)",
}


def draw_polar_polynomial_exactness(fig, n_mu, rows):
    k = np.array([r[0] for r in rows])
    err = np.array([r[3] for r in rows])
    within = np.array([r[4] for r in rows])
    err_plot = np.maximum(err, 1e-17)

    ax = fig.add_subplot(111)
    ax.semilogy(k[within], err_plot[within], "o-", color="tab:blue",
                label="within guaranteed range")
    ax.semilogy(k[~within], err_plot[~within], "o-", color="tab:red",
                label="beyond guaranteed range")
    ax.axvline(2 * n_mu - 1, color="gray", linestyle="--", linewidth=1,
               label=f"degree limit = 2*n_mu-1 = {2 * n_mu - 1}")
    ax.set_xlabel("degree k of mu^k")
    ax.set_ylabel("absolute error |numeric - exact|")
    ax.set_title(f"Polynomial exactness in mu (Gauss-Legendre, n_mu={n_mu})")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)


def draw_spherical_harmonic_orthonormality(fig, n_mu, n_phi, azimuthal, n_max, err_records):
    errs = np.array([r[4] for r in err_records])
    diag = np.array([r[5] for r in err_records])
    err_plot = np.maximum(errs, 1e-17)
    order = np.argsort(-err_plot)
    idx = np.arange(len(err_records))

    ax = fig.add_subplot(111)
    ax.semilogy(idx[~diag[order]], err_plot[order][~diag[order]], "o",
                color="tab:orange", markersize=3, label="off-diagonal pairs (n1,m1)!=(n2,m2)")
    ax.semilogy(idx[diag[order]], err_plot[order][diag[order]], "o",
                color="tab:blue", markersize=3, label="diagonal pairs (n1,m1)=(n2,m2)")
    ax.axhline(1e-10, color="gray", linestyle="--", linewidth=1,
               label="check threshold 1e-10")
    ax.set_xlabel(f"pairs (n1,m1,n2,m2), n1,n2<={n_max} -- sorted by decreasing error")
    ax.set_ylabel("absolute error")
    ax.set_title(f"Spherical harmonic orthonormality\n"
                 f"n_mu={n_mu}, n_phi={n_phi}, azimuthal='{azimuthal}'")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)


def draw_orthonormality_beyond_resolution(fig, n_mu, n_phi, curve, n_too_high, err):
    n = np.array([c[0] for c in curve])
    err_n = np.maximum(np.array([c[1] for c in curve]), 1e-17)

    ax = fig.add_subplot(111)
    ax.semilogy(n, err_n, "o-", color="tab:blue")
    ax.axvline(n_mu, color="gray", linestyle="--", linewidth=1,
               label=f"n_mu={n_mu} (expected resolution limit)")
    ax.axhline(1e-6, color="tab:red", linestyle=":", linewidth=1,
               label="'appreciable error' threshold 1e-6")
    ax.scatter([n_too_high], [max(err, 1e-17)], color="tab:red", zorder=5,
               label=f"n tested in the check = {n_too_high}")
    ax.set_xlabel("n (degree of Y_n^n tested against itself)")
    ax.set_ylabel("absolute error")
    ax.set_title(f"Quadrature resolution limit (n_mu={n_mu}, n_phi={n_phi})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)


def draw_quadrature_nodes_3d(fig, n_mu, n_phi, azimuthal, x, y, z, color_values, color_label):
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="lightgray", alpha=0.15, linewidth=0,
                     rstride=1, cstride=1, shade=False)

    if np.any(color_values < 0):
        vmax = np.max(np.abs(color_values))
        cmap, vmin, vmax_plot = "RdBu_r", -vmax, vmax
    else:
        cmap, vmin, vmax_plot = "viridis", None, None

    sc = ax.scatter(x, y, z, c=color_values, cmap=cmap, vmin=vmin, vmax=vmax_plot,
                     s=25, depthshade=True)
    fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, label=color_label)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z  (= mu)")
    ax.set_title(f"Product quadrature nodes on the unit sphere\n"
                 f"n_mu={n_mu}, n_phi={n_phi}, azimuthal='{azimuthal}' "
                 f"({n_mu * n_phi} nodes)")
    ax.set_box_aspect([1, 1, 1])


def draw_gauss_chebyshev_trapezoidal_equivalence(fig, n_phi, data):
    phi_trap = data["phi_trap"]
    phi_cheb = data["phi_cheb"]
    k_rows = data["k_rows"]
    k = np.array([r[0] for r in k_rows])
    err_trap = np.maximum(np.array([r[3] for r in k_rows]), 1e-17)
    err_cheb = np.maximum(np.array([r[4] for r in k_rows]), 1e-17)
    within = np.array([r[5] for r in k_rows])

    ax_nodes = fig.add_subplot(1, 2, 1, projection="polar")
    ax_err = fig.add_subplot(1, 2, 2)

    ax_nodes.scatter(phi_trap, np.ones_like(phi_trap), color="tab:blue", s=40,
                      label="trapezoidal nodes")
    ax_nodes.scatter(phi_cheb, np.ones_like(phi_cheb), color="tab:orange", s=40,
                      marker="x", label="Gauss-Chebyshev nodes")
    ax_nodes.set_yticklabels([])
    ax_nodes.set_title(f"Azimuthal nodes (n_phi={n_phi})", pad=20)
    ax_nodes.legend(loc="upper right", bbox_to_anchor=(1.3, 1.15), fontsize=8)

    ax_err.semilogy(k[within], err_trap[within], "o-", color="tab:blue", label="trapezoidal")
    ax_err.semilogy(k[within], err_cheb[within], "x--", color="tab:orange", label="Gauss-Chebyshev")
    ax_err.semilogy(k[~within], err_trap[~within], "o-", color="tab:blue", alpha=0.4)
    ax_err.semilogy(k[~within], err_cheb[~within], "x--", color="tab:orange", alpha=0.4)
    ax_err.axvline(n_phi - 1, color="gray", linestyle="--", linewidth=1,
                    label=f"degree limit = n_phi-1 = {n_phi - 1}")
    ax_err.set_xlabel("k in e^(i k phi)")
    ax_err.set_ylabel("absolute error")
    ax_err.set_title("Error vs degree k\n(solid = within range, faded = aliasing expected)")
    ax_err.legend(fontsize=8)
    ax_err.grid(True, which="both", alpha=0.3)

    fig.suptitle("Gauss-Chebyshev <-> trapezoidal equivalence")


class QuadratureCheckApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Sphere Quadrature Checks \u2014 Interactive GUI")
        self.root.geometry("1220x800")
        self.root.minsize(940, 620)

        self.field_vars: dict[str, tk.Variable] = {}
        self._last_check = None
        self._last_save = None  # no-arg(out_dir) -> path, set by each _run_xxx

        self._build_layout()
        self._on_check_changed()
        self._run_check()  # populate everything with sensible defaults on startup

    
    def _build_layout(self):
        bar = ttk.Frame(self.root, padding=8)
        bar.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(bar, text="Check:").pack(side=tk.LEFT)
        self.check_combo = ttk.Combobox(
            bar, state="readonly", width=34,
            values=[CHECK_LABELS[c] for c in _VALID_CHECKS],
        )
        self.check_combo.current(0)
        self.check_combo.pack(side=tk.LEFT, padx=(4, 12))
        self.check_combo.bind("<<ComboboxSelected>>", lambda e: self._on_check_changed())

        self.run_button = ttk.Button(bar, text="Run check", command=self._run_check)
        self.run_button.pack(side=tk.LEFT)

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)

        ttk.Label(bar, text="Output folder:").pack(side=tk.LEFT)
        self.out_dir_var = tk.StringVar(value="out")
        ttk.Entry(bar, textvariable=self.out_dir_var, width=22).pack(side=tk.LEFT, padx=(4, 4))
        ttk.Button(bar, text="Browse...", command=self._browse_out_dir).pack(side=tk.LEFT)
        ttk.Button(bar, text="Save results...", command=self._save_results).pack(side=tk.LEFT, padx=(12, 0))

        body = ttk.Frame(self.root)
        body.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        self.param_frame = ttk.LabelFrame(body, text="Parameters (= a section of input.txt)", padding=10)
        self.param_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 8))
        self.param_frame.columnconfigure(1, weight=1)

        right = ttk.Frame(body)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.notebook = ttk.Notebook(right)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.plot_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.plot_frame, text="Plot")
        self.fig = Figure(figsize=(8, 5.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        toolbar_frame = ttk.Frame(self.plot_frame)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

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

    def _current_check(self):
        return _VALID_CHECKS[self.check_combo.current()]

    def _add_label(self, row, text):
        ttk.Label(self.param_frame, text=text).grid(row=row, column=0, sticky="w", pady=(4, 0))

    def _add_int(self, row, name, label, default):
        self._add_label(row, label)
        var = tk.StringVar(value=str(default))
        ttk.Entry(self.param_frame, textvariable=var, width=14).grid(
            row=row, column=1, sticky="ew", pady=(4, 0), padx=(6, 0))
        self.field_vars[name] = var
        return row + 1

    def _add_choice(self, row, name, label, choices, default, on_change=None):
        self._add_label(row, label)
        var = tk.StringVar(value=default)
        combo = ttk.Combobox(self.param_frame, textvariable=var, values=list(choices),
                              state="readonly", width=18)
        combo.grid(row=row, column=1, sticky="ew", pady=(4, 0), padx=(6, 0))
        if on_change is not None:
            combo.bind("<<ComboboxSelected>>", lambda e: on_change())
        self.field_vars[name] = var
        return row + 1

    def _on_check_changed(self):
        prior = {name: var.get() for name, var in self.field_vars.items()}

        def default_for(name, fallback):
            return prior.get(name, fallback)

        for w in self.param_frame.winfo_children():
            w.destroy()
        self.field_vars = {}

        check = self._current_check()
        row = 0
        if check == "polar_polynomial_exactness":
            row = self._add_int(row, "n_mu", "n_mu (polar quadrature order)", default_for("n_mu", 6))
            row = self._add_int(row, "k_max", "k_max (blank behavior: 2*n_mu+2)", default_for("k_max", ""))

        elif check == "spherical_harmonic_orthonormality":
            row = self._add_int(row, "n_mu", "n_mu", default_for("n_mu", 6))
            row = self._add_int(row, "n_phi", "n_phi", default_for("n_phi", 13))
            row = self._add_int(row, "n_max", "n_max (max harmonic degree)", default_for("n_max", 4))
            row = self._add_choice(row, "azimuthal", "azimuthal rule", _VALID_AZIMUTHAL,
                                    default_for("azimuthal", _VALID_AZIMUTHAL[0]))

        elif check == "orthonormality_beyond_resolution":
            row = self._add_int(row, "n_mu", "n_mu", default_for("n_mu", 6))
            row = self._add_int(row, "n_phi", "n_phi", default_for("n_phi", 25))
            row = self._add_choice(row, "azimuthal", "azimuthal rule", _VALID_AZIMUTHAL,
                                    default_for("azimuthal", _VALID_AZIMUTHAL[0]))

        elif check == "gauss_chebyshev_trapezoidal_equivalence":
            row = self._add_int(row, "n_phi", "n_phi", default_for("n_phi", 13))

        else:  # quadrature_nodes_3d
            row = self._add_int(row, "n_mu", "n_mu", default_for("n_mu", 6))
            row = self._add_int(row, "n_phi", "n_phi", default_for("n_phi", 13))
            row = self._add_choice(row, "azimuthal", "azimuthal rule", _VALID_AZIMUTHAL,
                                    default_for("azimuthal", _VALID_AZIMUTHAL[0]))
            row = self._add_choice(row, "color_by", "color nodes by", _VALID_COLOR_BY,
                                    default_for("color_by", _VALID_COLOR_BY[0]),
                                    on_change=self._on_check_changed)
            if self.field_vars["color_by"].get() == "function":
                row = self._add_choice(row, "test_function", "test function", _VALID_TEST_FUNCTION,
                                        default_for("test_function", _VALID_TEST_FUNCTION[0]),
                                        on_change=self._on_check_changed)
                tf = self.field_vars["test_function"].get()
                if tf == "mu_power":
                    row = self._add_int(row, "k", "k (power of mu)", default_for("k", 4))
                elif tf == "spherical_harmonic":
                    row = self._add_int(row, "n", "n (harmonic degree)", default_for("n", 5))
                    row = self._add_int(row, "m", "m (harmonic order, |m|<=n)", default_for("m", 5))

    
    def _collect_int(self, name):
        raw = self.field_vars[name].get().strip()
        try:
            return int(raw)
        except ValueError:
            raise ValueError(f"'{name}' must be an integer (found '{raw}')")

    def _collect_optional_int(self, name):
        raw = self.field_vars[name].get().strip()
        if not raw:
            return None
        try:
            return int(raw)
        except ValueError:
            raise ValueError(f"'{name}' must be an integer or blank (found '{raw}')")

    def _require_positive(self, name, value):
        if value <= 0:
            raise ValueError(f"'{name}' must be a positive integer")

    def _run_check(self):
        check = self._current_check()
        try:
            if check == "polar_polynomial_exactness":
                n_mu = self._collect_int("n_mu")
                self._require_positive("n_mu", n_mu)
                k_max = self._collect_optional_int("k_max")
                self._exec(check, lambda: self._run_polar_polynomial(n_mu, k_max))

            elif check == "spherical_harmonic_orthonormality":
                n_mu = self._collect_int("n_mu")
                self._require_positive("n_mu", n_mu)
                n_phi = self._collect_int("n_phi")
                self._require_positive("n_phi", n_phi)
                n_max = self._collect_int("n_max")
                if n_max < 0:
                    raise ValueError("'n_max' must be >= 0")
                azimuthal = self.field_vars["azimuthal"].get()
                self._exec(check, lambda: self._run_spherical_orthonormality(n_mu, n_phi, n_max, azimuthal))

            elif check == "orthonormality_beyond_resolution":
                n_mu = self._collect_int("n_mu")
                self._require_positive("n_mu", n_mu)
                n_phi = self._collect_int("n_phi")
                self._require_positive("n_phi", n_phi)
                azimuthal = self.field_vars["azimuthal"].get()
                self._exec(check, lambda: self._run_beyond_resolution(n_mu, n_phi, azimuthal))

            elif check == "gauss_chebyshev_trapezoidal_equivalence":
                n_phi = self._collect_int("n_phi")
                self._require_positive("n_phi", n_phi)
                self._exec(check, lambda: self._run_gc_trap_equivalence(n_phi))

            else:  # quadrature_nodes_3d
                n_mu = self._collect_int("n_mu")
                self._require_positive("n_mu", n_mu)
                n_phi = self._collect_int("n_phi")
                self._require_positive("n_phi", n_phi)
                azimuthal = self.field_vars["azimuthal"].get()
                color_by = self.field_vars["color_by"].get()
                test_function = k_val = n_val = m_val = None
                if color_by == "function":
                    test_function = self.field_vars["test_function"].get()
                    if test_function == "mu_power":
                        k_val = self._collect_int("k")
                    elif test_function == "spherical_harmonic":
                        n_val = self._collect_int("n")
                        m_val = self._collect_int("m")
                        if abs(m_val) > n_val:
                            raise ValueError("'m' must satisfy |m| <= n")
                self._exec(check, lambda: self._run_quadrature_nodes_3d(
                    n_mu, n_phi, azimuthal, color_by, test_function, k_val, n_val, m_val))

        except ValueError as e:
            messagebox.showerror("Invalid input", str(e))

    def _exec(self, check, body_fn):
        self.run_button.config(state=tk.DISABLED)
        self.status_var.set(f"Running {check}...")
        self.root.update_idletasks()

        buf = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buf
        try:
            self.fig.clear()
            body_fn()
            self.canvas.draw()
        except Exception as e:  # noqa: BLE001 -- surface any error to the user
            sys.stdout = old_stdout
            self.run_button.config(state=tk.NORMAL)
            messagebox.showerror("Check failed", str(e))
            self.status_var.set("Error.")
            return
        finally:
            sys.stdout = old_stdout

        self._set_report_text(buf.getvalue() or "(no textual report for this check)")
        self._last_check = check
        self.run_button.config(state=tk.NORMAL)
        self.status_var.set(f"Done: {check}")
        self.notebook.select(self.plot_frame)

    def _set_report_text(self, text):
        self.report_text.config(state="normal")
        self.report_text.delete("1.0", tk.END)
        self.report_text.insert(tk.END, text)
        self.report_text.config(state="disabled")

    def _run_polar_polynomial(self, n_mu, k_max):
        rows, max_err_within = checks.polar_polynomial_errors(n_mu=n_mu, k_max=k_max)
        report.report_polar_polynomial_exactness(n_mu=n_mu, rows=rows, max_err_within=max_err_within)
        if max_err_within >= 1e-12:
            print("\n[WARNING] expected exact to machine precision (max_err_within < 1e-12)!")
        draw_polar_polynomial_exactness(self.fig, n_mu, rows)
        self._last_save = lambda out_dir: self._save(
            out_dir, "polar_polynomial_exactness",
            lambda path: plots.plot_polar_polynomial_exactness(n_mu=n_mu, rows=rows, out_path=path))

    def _run_spherical_orthonormality(self, n_mu, n_phi, n_max, azimuthal):
        worst_err, n_pairs, err_records = checks.spherical_harmonic_orthonormality_errors(
            n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal, n_max=n_max)
        report.report_spherical_harmonic_orthonormality(
            n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal, n_max=n_max,
            worst_err=worst_err, n_pairs_checked=n_pairs)
        if worst_err >= 1e-10:
            print("\n[WARNING] expected an error at machine precision within the resolution!")
        draw_spherical_harmonic_orthonormality(self.fig, n_mu, n_phi, azimuthal, n_max, err_records)
        self._last_save = lambda out_dir: self._save(
            out_dir, "spherical_harmonic_orthonormality",
            lambda path: plots.plot_spherical_harmonic_orthonormality(
                n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal, n_max=n_max,
                err_records=err_records, out_path=path))

    def _run_beyond_resolution(self, n_mu, n_phi, azimuthal):
        n_too_high, err_beyond, curve = checks.orthonormality_beyond_resolution_error(
            n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal)
        report.report_orthonormality_beyond_resolution(
            n_mu=n_mu, n_phi=n_phi, n_too_high=n_too_high, err=err_beyond)
        if err_beyond <= 1e-6:
            print("\n[WARNING] expected an appreciable error beyond the resolution!")
        draw_orthonormality_beyond_resolution(self.fig, n_mu, n_phi, curve, n_too_high, err_beyond)
        self._last_save = lambda out_dir: self._save(
            out_dir, "orthonormality_beyond_resolution",
            lambda path: plots.plot_orthonormality_beyond_resolution(
                n_mu=n_mu, n_phi=n_phi, curve=curve, n_too_high=n_too_high,
                err=err_beyond, out_path=path))

    def _run_gc_trap_equivalence(self, n_phi):
        data = checks.gauss_chebyshev_trapezoidal_equivalence(n_phi=n_phi)
        report.report_gauss_chebyshev_trapezoidal_equivalence(n_phi=n_phi, data=data)
        draw_gauss_chebyshev_trapezoidal_equivalence(self.fig, n_phi, data)
        self._last_save = lambda out_dir: self._save(
            out_dir, "gauss_chebyshev_trapezoidal_equivalence",
            lambda path: plots.plot_gauss_chebyshev_trapezoidal_equivalence(
                n_phi=n_phi, data=data, out_path=path))

    def _run_quadrature_nodes_3d(self, n_mu, n_phi, azimuthal, color_by, test_function, k, n, m):
        quad, x, y, z, color_values, color_label = checks.quadrature_node_view(
            n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal, color_by=color_by,
            test_function=test_function, k=k, n=n, m=m)
        report.report_quadrature_nodes_3d(
            n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal, quad=quad,
            color_by=color_by, color_label=color_label)
        draw_quadrature_nodes_3d(self.fig, n_mu, n_phi, azimuthal, x, y, z, color_values, color_label)
        self._last_save = lambda out_dir: self._save(
            out_dir, "quadrature_nodes_3d",
            lambda path: plots.plot_quadrature_nodes_3d(
                n_mu=n_mu, n_phi=n_phi, azimuthal=azimuthal,
                x=x, y=y, z=z, color_values=color_values, color_label=color_label, out_path=path))

    def _save(self, out_dir, check_name, plot_fn):
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, f"plot_{check_name}.png")
        plot_fn(path)
        return path

    def _browse_out_dir(self):
        directory = filedialog.askdirectory(title="Choose output folder")
        if directory:
            self.out_dir_var.set(directory)

    def _save_results(self):
        if self._last_save is None:
            messagebox.showinfo("Nothing to save", "Run a check first.")
            return
        out_dir = self.out_dir_var.get().strip() or "out"
        try:
            plot_path = self._last_save(out_dir)
            log_path = os.path.join(out_dir, f"report_{self._last_check}.txt")
            with open(log_path, "w", encoding="utf-8") as f:
                f.write(self.report_text.get("1.0", tk.END))
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Save failed", str(e))
            return
        self.status_var.set(f"Saved report and plot to {os.path.abspath(out_dir)}")
        messagebox.showinfo("Saved", f"Wrote:\n{plot_path}\n{log_path}")


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    QuadratureCheckApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
