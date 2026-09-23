"""
gui/app.py
----------
Tkinter desktop GUI for the milne_pn package: interactive plot and data
visualization for the PN solution of the Milne problem.

Run it with:
    python -m milne_pn.gui
    # or, after `pip install -e .`:
    milne-gui

Layout
------
  * Top control bar: pick a PN order and solve it; run a convergence
    sweep up to a chosen N_max (runs in a background thread so the UI
    stays responsive); export CSVs / save the current plot.
  * A tabbed notebook with four embedded matplotlib figures (Flux
    Profile, Convergence, Gibbs Ringing, Runge Phenomenon) and a Data
    Table tab (a sortable ttk.Treeview) showing the convergence results
    numerically.
  * A status bar reporting the latest z0 and its error against Case's
    exact value.

This module only imports the plain, headless algorithms/benchmarks
modules -- it holds all the Tkinter/matplotlib-embedding logic itself, so
the rest of the package stays GUI-independent.
"""

from __future__ import annotations

import os
import queue
import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from ..algorithms.milne_solver import MilneSolver
from ..algorithms.runge_analysis import (
    reconstruct_boundary_flux,
    gibbs_overshoot,
    runge_interpolation_demo,
)
from ..benchmarks.benchmark import convergence_benchmark
from ..utils.config import load_config
from ..utils.io_utils import write_csv

Z0_EXACT = 0.710446
DEFAULT_GIBBS_ORDERS = [3, 9, 21, 41]
DEFAULT_RUNGE_DEGREES = [5, 10, 15, 20]


# ---------------------------------------------------------------------------
# small reusable widget: one notebook tab holding an embedded matplotlib Figure
# ---------------------------------------------------------------------------
class PlotTab:
    def __init__(self, notebook, title, figsize=(9, 5)):
        self.frame = ttk.Frame(notebook)
        notebook.add(self.frame, text=title)

        self.fig = Figure(figsize=figsize, dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)

        toolbar_frame = ttk.Frame(self.frame)
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self.toolbar.update()

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def clear(self):
        self.fig.clear()

    def draw(self):
        self.fig.tight_layout()
        self.canvas.draw()


# ---------------------------------------------------------------------------
# main application
# ---------------------------------------------------------------------------
class MilneApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Milne PN Solver \u2014 Plot & Data Visualization")
        self.root.geometry("1150x780")
        self.root.minsize(820, 560)

        self.cfg = load_config()
        self._task_queue: "queue.Queue" = queue.Queue()

        self._last_solver: MilneSolver | None = None
        self._convergence_result: dict | None = None
        self._gibbs_rows: list | None = None
        self._gibbs_orders = list(self.cfg.get("runge_analysis", {}).get("orders", DEFAULT_GIBBS_ORDERS))

        self._build_menu()
        self._build_control_bar()
        self._build_tabs()
        self._build_status_bar()

        self._poll_job = self.root.after(100, self._poll_queue)
        self.root.bind("<Destroy>", self._on_destroy, add="+")

        # populate everything with sensible defaults on startup
        self._solve_single()
        self._plot_runge(DEFAULT_RUNGE_DEGREES)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_menu(self):
        menubar = tk.Menu(self.root)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Export CSVs...", command=self._export_csvs)
        file_menu.add_command(label="Save Current Plot...", command=self._save_current_plot)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.root.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

    def _build_control_bar(self):
        bar = ttk.Frame(self.root, padding=8)
        bar.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(bar, text="PN order N (odd):").pack(side=tk.LEFT)
        self.n_var = tk.IntVar(value=int(self.cfg["solver"]["n_plot"]))
        self.n_spin = ttk.Spinbox(bar, from_=1, to=299, increment=2, textvariable=self.n_var, width=6,
                                   command=self._solve_single)
        self.n_spin.pack(side=tk.LEFT, padx=(4, 4))
        self.n_spin.bind("<Return>", lambda e: self._solve_single())

        self.solve_button = ttk.Button(bar, text="Solve", command=self._solve_single)
        self.solve_button.pack(side=tk.LEFT, padx=(4, 0))

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)

        ttk.Label(bar, text="Convergence N max:").pack(side=tk.LEFT)
        self.nmax_var = tk.IntVar(value=int(self.cfg["solver"]["n_max"]))
        self.nmax_spin = ttk.Spinbox(bar, from_=3, to=299, increment=2, textvariable=self.nmax_var, width=6)
        self.nmax_spin.pack(side=tk.LEFT, padx=(4, 4))

        self.conv_button = ttk.Button(bar, text="Run Convergence Study", command=self._run_convergence)
        self.conv_button.pack(side=tk.LEFT, padx=(4, 0))

        self.progress = ttk.Progressbar(bar, mode="indeterminate", length=140)
        self.progress.pack(side=tk.LEFT, padx=12)

        ttk.Separator(bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)

        ttk.Button(bar, text="Export CSVs...", command=self._export_csvs).pack(side=tk.LEFT)
        ttk.Button(bar, text="Save Current Plot...", command=self._save_current_plot).pack(side=tk.LEFT, padx=(6, 0))

    def _build_tabs(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 8))

        self.tab_flux = PlotTab(self.notebook, "Flux Profile")
        self.tab_conv = PlotTab(self.notebook, "Convergence")
        self.tab_gibbs = PlotTab(self.notebook, "Gibbs Ringing")
        self.tab_runge = PlotTab(self.notebook, "Runge Phenomenon")
        self._plot_tabs = [self.tab_flux, self.tab_conv, self.tab_gibbs, self.tab_runge]

        self._build_data_tab()

    def _build_data_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Data Table")
        self._data_tab_index = len(self._plot_tabs)  # index within notebook children

        columns = ("N", "z0", "error", "modes", "max_imag")
        headers = {"N": "N", "z0": "z0", "error": "error", "modes": "# modes", "max_imag": "max|Im(\u03c9)|"}
        widths = {"N": 70, "z0": 150, "error": 150, "modes": 90, "max_imag": 130}

        self.tree = ttk.Treeview(frame, columns=columns, show="headings", height=22)
        for c in columns:
            self.tree.heading(c, text=headers[c], command=lambda col=c: self._sort_tree(col))
            self.tree.column(c, width=widths[c], anchor=tk.E)

        vsb = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4, 0), pady=4)
        vsb.pack(side=tk.RIGHT, fill=tk.Y, pady=4)

    def _build_status_bar(self):
        self.status_var = tk.StringVar(value="Ready.")
        bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W, padding=(6, 2))
        bar.pack(side=tk.BOTTOM, fill=tk.X)

    # ------------------------------------------------------------------
    # actions
    # ------------------------------------------------------------------
    def _normalize_odd(self, var: tk.IntVar) -> int:
        try:
            v = int(var.get())
        except (tk.TclError, ValueError):
            v = 1
        if v < 1:
            v = 1
        if v % 2 == 0:
            v += 1
        var.set(v)
        return v

    def _solve_single(self):
        N = self._normalize_odd(self.n_var)
        try:
            solver = MilneSolver(N)
        except Exception as e:  # noqa: BLE001 - surface any solver error to the user
            messagebox.showerror("Solve failed", str(e))
            return

        self._last_solver = solver
        self._plot_flux(solver)

        # keep the Gibbs-ringing tab's order list current (always include N)
        orders = sorted(set(self._gibbs_orders) | {N})
        self._plot_gibbs(orders)

        err = solver.z0 - Z0_EXACT
        self.status_var.set(
            f"N={N}   z0={solver.z0:.8f}   error={err:+.2e}   "
            f"discrete modes={solver.num_discrete_modes}   "
            f"max|Im(\u03c9)|={solver.max_imag_residual:.1e}"
        )

    def _plot_flux(self, solver: MilneSolver):
        tab = self.tab_flux
        tab.clear()
        ax = tab.fig.add_subplot(111)

        x = np.linspace(0, 6, 300)
        phi0 = solver.phi0(x)
        ax.plot(x, phi0, color="steelblue", lw=2, label=rf"$\phi_0(x)$ ($P_{{{solver.N}}}$)")
        ax.plot(x, x + solver.z0, color="crimson", ls="--", lw=1.5, label=r"asymptote $x+z_0$")

        xline = np.linspace(-solver.z0, 0, 50)
        ax.plot(xline, xline + solver.z0, color="crimson", ls=":", lw=1.5)
        ax.axvline(0, color="gray", lw=0.8)
        ax.axhline(0, color="gray", lw=0.8)
        ax.scatter([-solver.z0], [0], color="crimson", zorder=5)
        ax.annotate(rf"$-z_0={-solver.z0:.4f}$", (-solver.z0, 0),
                    textcoords="offset points", xytext=(8, -14))

        ax.set_xlim(-1, 6)
        ax.set_xlabel("x  [mean free paths]")
        ax.set_ylabel(r"$\phi_0(x)$")
        ax.set_title(f"Milne problem scalar flux, $P_{{{solver.N}}}$ approximation")
        ax.legend()
        ax.grid(alpha=0.3)
        tab.draw()

    def _plot_gibbs(self, orders):
        tab = self.tab_gibbs
        tab.clear()

        ax1 = tab.fig.add_subplot(121)
        colors = matplotlib.colormaps["viridis"](np.linspace(0.15, 0.85, len(orders)))
        rows = []
        for N, color in zip(orders, colors):
            solver = MilneSolver(N)
            mu, psi0 = reconstruct_boundary_flux(solver, n_mu=1000)
            ax1.plot(mu, psi0, color=color, lw=1.3, label=rf"$P_{{{N}}}$")
            rows.append(gibbs_overshoot(solver, n_mu=1000))
        self._gibbs_rows = rows

        ax1.axvline(0, color="black", lw=0.8)
        ax1.axhline(0, color="black", lw=0.6)
        ax1.axvspan(0, 1, color="red", alpha=0.05)
        ax1.set_xlabel(r"$\mu$")
        ax1.set_ylabel(r"$\psi(0,\mu)$")
        ax1.set_title("Gibbs ringing at $\\mu=0$\n(shaded: should be exactly 0)")
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)

        ax2 = tab.fig.add_subplot(122)
        rel = [r["relative_overshoot"] * 100 for r in rows]
        ax2.plot(orders, rel, "o-", color="firebrick")
        ax2.set_xlabel("PN order N")
        ax2.set_ylabel("overshoot in $\\mu>0$ [%]")
        ax2.set_title("Overshoot plateau\n(classic Gibbs signature)")
        ax2.grid(alpha=0.3)
        ax2.set_ylim(bottom=0)

        tab.draw()

    def _plot_runge(self, degrees):
        tab = self.tab_runge
        tab.clear()
        result = runge_interpolation_demo(degrees=degrees, n_eval=1000)
        x = result["x_eval"]

        ax1 = tab.fig.add_subplot(121)
        ax1.plot(x, result["y_true"], "k-", lw=2, label=r"$f(x)=1/(1{+}25x^2)$")
        colors = matplotlib.colormaps["viridis"](np.linspace(0.15, 0.85, len(degrees)))
        for d, color in zip(degrees, colors):
            ax1.plot(x, result["interpolants"][d], color=color, lw=1.2, label=f"degree {d}")
        ax1.set_ylim(-1, 2)
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_title("Classical Runge phenomenon\n(equispaced-node interpolation)")
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)

        ax2 = tab.fig.add_subplot(122)
        degs = sorted(result["max_abs_error"].keys())
        errs = [result["max_abs_error"][d] for d in degs]
        ax2.semilogy(degs, errs, "o-", color="firebrick")
        ax2.set_xlabel("interpolation degree")
        ax2.set_ylabel("max abs. error")
        ax2.set_title("Error GROWS with degree")
        ax2.grid(alpha=0.3, which="both")

        tab.draw()

    def _run_convergence(self):
        nmax = self._normalize_odd(self.nmax_var)
        self.conv_button.config(state=tk.DISABLED)
        self.solve_button.config(state=tk.DISABLED)
        self.progress.start(10)
        self.status_var.set(f"Running convergence study up to N={nmax}...")

        thread = threading.Thread(target=self._convergence_worker, args=(nmax,), daemon=True)
        thread.start()

    def _convergence_worker(self, nmax):
        try:
            result = convergence_benchmark(n_max=nmax, z0_exact=Z0_EXACT)
            self._task_queue.put(("convergence_done", result))
        except Exception as e:  # noqa: BLE001
            self._task_queue.put(("error", str(e)))

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self._task_queue.get_nowait()
                if kind == "convergence_done":
                    self._on_convergence_done(payload)
                elif kind == "error":
                    self.progress.stop()
                    self.conv_button.config(state=tk.NORMAL)
                    self.solve_button.config(state=tk.NORMAL)
                    messagebox.showerror("Error", payload)
        except queue.Empty:
            pass
        self._poll_job = self.root.after(100, self._poll_queue)

    def _on_destroy(self, event):
        # avoid a stray after()-callback firing against an already-destroyed
        # root (harmless in normal single-instance use, but tidy this up so
        # embedding/testing code that creates+destroys multiple app
        # instances in one process doesn't see spurious Tcl errors)
        if event.widget is self.root and getattr(self, "_poll_job", None) is not None:
            try:
                self.root.after_cancel(self._poll_job)
            except tk.TclError:
                pass
            self._poll_job = None

    def _on_convergence_done(self, result):
        self.progress.stop()
        self.conv_button.config(state=tk.NORMAL)
        self.solve_button.config(state=tk.NORMAL)
        self._convergence_result = result

        self._plot_convergence(result)
        self._fill_data_table(result)
        self.notebook.select(self.tab_conv.frame)

        self.status_var.set(
            f"Convergence study done: N up to {result['N'][-1]}, "
            f"z0={result['z0'][-1]:.8f} (error {result['error'][-1]:+.2e})"
        )

    def _plot_convergence(self, result):
        tab = self.tab_conv
        tab.clear()

        ax1 = tab.fig.add_subplot(121)
        ax1.axhline(Z0_EXACT, color="crimson", ls="--", lw=1.3, label=f"Case exact: {Z0_EXACT}")
        ax1.plot(result["N"], result["z0"], "o-", color="steelblue", ms=4, label=r"$P_N$ result")
        ax1.set_xlabel("PN order N")
        ax1.set_ylabel(r"$z_0$ [mfp]")
        ax1.set_title(r"Convergence of $z_0$")
        ax1.legend()
        ax1.grid(alpha=0.3)

        ax2 = tab.fig.add_subplot(122)
        err = np.abs(result["error"])
        ax2.semilogy(result["N"], np.maximum(err, 1e-16), "o-", color="darkorange", ms=4)
        ax2.set_xlabel("PN order N")
        ax2.set_ylabel(r"$|z_0^{(N)}-z_0^{exact}|$")
        ax2.set_title("Convergence error (log scale)")
        ax2.grid(alpha=0.3, which="both")

        tab.draw()

    def _fill_data_table(self, result):
        for row in self.tree.get_children():
            self.tree.delete(row)
        for N, z0, err, modes, mi in zip(
            result["N"], result["z0"], result["error"],
            result["num_discrete_modes"], result["max_imag_residual"],
        ):
            self.tree.insert("", tk.END, values=(
                int(N), f"{z0:.8f}", f"{err:+.2e}", int(modes), f"{mi:.1e}",
            ))

    def _sort_tree(self, col):
        items = [(self.tree.set(k, col), k) for k in self.tree.get_children("")]
        try:
            items.sort(key=lambda t: float(t[0].replace("+", "")))
        except ValueError:
            items.sort()
        for index, (_, k) in enumerate(items):
            self.tree.move(k, "", index)

    # ------------------------------------------------------------------
    # export / save
    # ------------------------------------------------------------------
    def _export_csvs(self):
        if self._convergence_result is None and self._last_solver is None:
            messagebox.showinfo("Nothing to export", "Solve a PN order or run a convergence study first.")
            return

        directory = filedialog.askdirectory(title="Choose export folder")
        if not directory:
            return

        written = []
        try:
            if self._convergence_result is not None:
                r = self._convergence_result
                path = os.path.join(directory, "convergence.csv")
                write_csv(path, {
                    "N": r["N"], "z0": r["z0"], "error": r["error"],
                    "num_discrete_modes": r["num_discrete_modes"],
                    "max_imag_residual": r["max_imag_residual"],
                })
                written.append(path)

            if self._last_solver is not None:
                solver = self._last_solver
                x = np.linspace(0, 6, 401)
                path = os.path.join(directory, "flux_profile.csv")
                write_csv(path, {"x": x, "phi0": solver.phi0(x), "phi0_asymptotic": x + solver.z0})
                written.append(path)

            if self._gibbs_rows:
                path = os.path.join(directory, "gibbs_overshoot.csv")
                write_csv(path, {
                    "N": np.array([r["N"] for r in self._gibbs_rows]),
                    "max_overshoot_abs": np.array([r["max_overshoot_abs"] for r in self._gibbs_rows]),
                    "outgoing_scale": np.array([r["outgoing_scale"] for r in self._gibbs_rows]),
                    "relative_overshoot": np.array([r["relative_overshoot"] for r in self._gibbs_rows]),
                })
                written.append(path)
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Export failed", str(e))
            return

        self.status_var.set(f"Exported {len(written)} CSV file(s) to {directory}")
        messagebox.showinfo("Export complete", "Wrote:\n" + "\n".join(written))

    def _current_plot_tab(self):
        selected = self.notebook.select()
        for tab in self._plot_tabs:
            if str(tab.frame) == str(selected):
                return tab
        return None

    def _save_current_plot(self):
        tab = self._current_plot_tab()
        if tab is None:
            messagebox.showinfo("Nothing to save", "Select a plot tab first (Data Table has no figure).")
            return
        path = filedialog.asksaveasfilename(
            title="Save current plot as...", defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            tab.fig.savefig(path, dpi=150)
        except Exception as e:  # noqa: BLE001
            messagebox.showerror("Save failed", str(e))
            return
        self.status_var.set(f"Saved plot to {path}")

    def _show_about(self):
        messagebox.showinfo(
            "About",
            "Milne PN Solver \u2014 Plot & Data Visualization\n\n"
            "PN (Legendre) approximation to the Milne problem, with a "
            "Runge/Gibbs-phenomenon analysis of the PN convergence.\n\n"
            f"Case's exact extrapolation distance: z0 = {Z0_EXACT} mfp",
        )


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    MilneApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
