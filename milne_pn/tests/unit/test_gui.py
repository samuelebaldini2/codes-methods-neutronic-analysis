"""
Unit tests for milne_pn.gui.app.

These need a real (or virtual, e.g. Xvfb) X display to create a Tk root
window, which typical CI runners don't have by default. The whole module
is skipped cleanly if tkinter can't be imported or no Tk root can be
created, rather than failing the suite.
"""

import pytest

tk = pytest.importorskip("tkinter")

try:
    _root = tk.Tk()
    _root.destroy()
    _HAS_DISPLAY = True
except Exception:
    _HAS_DISPLAY = False

pytestmark = pytest.mark.skipif(not _HAS_DISPLAY, reason="no usable Tk display available")


@pytest.fixture
def app():
    from tkinter import ttk
    from milne_pn.gui.app import MilneApp

    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    application = MilneApp(root)
    root.update()
    yield application
    root.destroy()


def test_app_solves_on_startup(app):
    assert app._last_solver is not None
    assert abs(app._last_solver.z0 - 0.710446) < 0.01


def test_normalize_odd_rounds_up(app):
    app.n_var.set(20)
    assert app._normalize_odd(app.n_var) == 21


def test_normalize_odd_leaves_odd_values_alone(app):
    app.n_var.set(15)
    assert app._normalize_odd(app.n_var) == 15


def test_normalize_odd_floors_at_one(app):
    app.n_var.set(-5)
    assert app._normalize_odd(app.n_var) == 1


def test_solve_single_updates_status(app):
    app.n_var.set(9)
    app._solve_single()
    assert "N=9" in app.status_var.get()
    assert app._last_solver.N == 9


def test_all_plot_tabs_are_selectable(app):
    for tab in app._plot_tabs:
        app.notebook.select(tab.frame)
        assert app._current_plot_tab() is tab


def test_data_table_tab_has_no_plot(app):
    for i in range(len(app.notebook.tabs())):
        if app.notebook.tab(i, "text") == "Data Table":
            app.notebook.select(i)
            break
    assert app._current_plot_tab() is None


def test_convergence_study_populates_table(app):
    app.nmax_var.set(11)
    app._run_convergence()
    import time
    deadline = time.time() + 20
    while app._convergence_result is None and time.time() < deadline:
        app._poll_queue()
        app.root.update()
        time.sleep(0.05)
    assert app._convergence_result is not None
    assert len(app.tree.get_children()) == 6  # N = 1,3,5,7,9,11
