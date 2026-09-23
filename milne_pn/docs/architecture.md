# Architecture

## Directory layout

```text
milne-pn/
|-- src/milne_pn/
|   |-- __init__.py           # exposes MilneSolver, PNSystem at top level
|   |-- algorithms/            # the mathematical core
|   |   |-- pn_system.py        # PNSystem: builds A, D, H
|   |   |-- milne_solver.py     # MilneSolver: eigenmodes + Marshak closure -> z0
|   |   `-- runge_analysis.py   # Gibbs-ringing / Runge-phenomenon analysis
|   |-- benchmarks/
|   |   `-- benchmark.py        # convergence_benchmark(), timing_benchmark()
|   |-- visualization/
|   |   |-- convergence_visualizer.py
|   |   |-- flux_visualizer.py
|   |   `-- angular_flux_visualizer.py   # Gibbs ringing + Runge demo plots
|   |-- gui/
|   |   |-- app.py              # Tkinter desktop app (MilneApp, main())
|   |   `-- __main__.py         # `python -m milne_pn.gui`
|   `-- utils/
|       |-- quadrature.py       # Gauss-Legendre + Legendre polynomials (numpy)
|       |-- config.py           # config/config.yaml loader with defaults
|       `-- io_utils.py         # dependency-free CSV read/write
|-- tests/
|   |-- unit/                   # one test module per algorithms/utils/gui module
|   |-- integration/            # full pipeline, against Case exact z0
|   `-- performance/            # timing sanity checks
|-- examples/                   # runnable, self-contained scripts
|-- notebooks/                  # interactive walkthroughs (mirror the examples)
|-- scripts/                    # thin CLIs over the package, config-aware
|-- data/{input,output}/        # generated CSVs land in data/output/
|-- docs/                       # this file, algorithms.md, performance.md
|-- config/config.yaml
|-- .github/workflows/tests.yml
|-- pyproject.toml / requirements.txt
`-- README.md
```

## Design principles

Layered, one-directional dependencies. utils has no dependency on
anything else in the package. algorithms depends only on utils.
benchmarks and visualization depend on algorithms (and visualization
additionally on algorithms.runge_analysis for the Gibbs plots). gui
depends on algorithms, benchmarks, and visualization, and is the only
subpackage that imports tkinter. scripts, examples, and notebooks sit on
top and depend on all of the above. Nothing lower in this list ever
imports from something higher.

All numerics through numpy built-ins. No hand-rolled linear algebra:
np.linalg.solve, np.linalg.eig, and np.polynomial.legendre do the
matrix-solve, eigenproblem, and quadrature work respectively. (A separate,
dependency-free C++ port of this same method exists as a companion
project and does implement its own dense linear algebra from scratch,
since it has no access to a library like Eigen/LAPACK in its target
environment -- but that is a constraint specific to that port, not a
design choice replicated here.)

The defective omega=0 mode is never asked of the eigensolver. np.linalg.eig
cannot reliably return a Jordan chain for a non-diagonalizable (defective)
eigenvalue. Rather than fight the eigensolver, MilneSolver derives that
mode exact closed form analytically once (see docs/algorithms.md) and
uses np.linalg.eig only for the remaining, well-conditioned, genuinely
diagonalizable nonzero eigenvalues.

Visualization -- both the matplotlib figure functions AND the Tkinter GUI --
is a thin, optional layer. algorithms/ and benchmarks/ have no plotting
dependency; only visualization/ and gui/ import matplotlib, and only gui/
imports tkinter. This keeps `pip install milne-pn` (core) lightweight, with
`pip install milne-pn[gui]` (or [dev], which also pulls in pytest) opt-in
for anything that renders a figure or opens a window. The GUI reuses the
exact same headless functions as the CLI scripts and the plain matplotlib
visualizers (MilneSolver, convergence_benchmark, gibbs_overshoot,
runge_interpolation_demo) -- it adds no separate computational path, only
a Tkinter shell around the same library calls.

Config over hardcoding. PN orders, quadrature resolution, tolerances, and
output paths live in config/config.yaml, loaded by utils/config.py. Every
script (scripts/*.py) and the GUI itself read this file (for their
startup defaults) but accept overrides on top of it (CLI flags for
scripts, the N/N-max spinboxes for the GUI); direct package/API use
(examples/, notebooks/) can ignore config entirely.

CSV as the interchange format. benchmarks/benchmark.py, scripts/generate_data.py,
and the GUI CSV export action all write plain CSVs (utils/io_utils.py,
stdlib csv + numpy, no pandas dependency) into data/output/ (or wherever
the user picks in the GUI export dialog). visualization/ and the notebooks
consume these CSVs rather than recomputing, so a plot can always be
regenerated without re-running the solver, and the schema is stable
enough for an external tool (or another language port) to consume too.

## The GUI specifically

`gui/app.py` holds a `PlotTab` helper (one notebook tab = one embedded
`matplotlib.figure.Figure` + `FigureCanvasTkAgg` + `NavigationToolbar2Tk`)
and the `MilneApp` class, which wires up:

  * a single-N solve (`MilneSolver`, synchronous -- sub-second even at
    N~150, so no threading needed),
  * a convergence sweep (`convergence_benchmark`, run on a background
    `threading.Thread` with results handed back through a `queue.Queue`
    and drained by a `root.after()` poll loop, so the Tk event loop is
    never blocked),
  * the Gibbs-ringing and Runge-phenomenon plots (`gibbs_overshoot`,
    `reconstruct_boundary_flux`, `runge_interpolation_demo` -- the same
    functions the notebooks and `visualization/angular_flux_visualizer.py`
    use), and
  * a `ttk.Treeview` data table kept in sync with the latest convergence
    run, sortable by clicking any column header.

No GUI code re-implements any solver or plotting logic; it only calls into
`algorithms/`, `benchmarks/`, and reuses the same plotting recipes as
`visualization/` (adapted to draw onto an embedded `Figure` instead of
saving directly to a file).

## Data flow

```
config/config.yaml
        |
        v
scripts/generate_data.py  --->  algorithms/{pn_system,milne_solver,runge_analysis}.py
        |                              |                        ^
        v                              v                        |
data/output/*.csv            benchmarks/benchmark.py      gui/app.py (same calls,
        |                                                  rendered into Tk tabs
        v                                                  instead of files)
visualization/*.py  --->  PNG plots (examples/output/, notebooks, or ad hoc)
```

## Testing strategy

  * tests/unit/ exercises each module (PNSystem, MilneSolver, quadrature,
    and the GUI) in isolation. The GUI tests (tests/unit/test_gui.py) spin
    up and tear down a real (or Xvfb virtual) Tk root per test via a
    pytest fixture, and the whole module is skipped cleanly (not failed)
    if no usable display is available -- most CI runners do not have one
    by default.
  * tests/integration/ exercises the full pipeline end to end, including
    the Runge/Gibbs analysis, checking convergence toward Case 0.710446
    and confirming the Gibbs-overshoot plateau.
  * tests/performance/ guards against gross performance regressions
    rather than asserting tight timings, since CI runners vary too much
    for strict wall-clock assertions.
