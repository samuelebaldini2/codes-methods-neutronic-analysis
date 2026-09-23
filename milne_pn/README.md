# milne-pn

A **P_N (Legendre / spherical-harmonics) approximation** solver for the
**Milne problem**: a semi-infinite, purely-scattering half-space
(x ∈ [0, ∞), c = 1, no absorption, no internal source), with the source
effectively at x → +∞. Computes the **extrapolation distance z₀** and its
convergence, as the odd PN order N increases, toward Case's exact value
**z₀ = 0.7104** mean free paths — and explains *why* that convergence is
slow: a **Runge/Gibbs-phenomenon** analysis of the truncated angular-flux
reconstruction, built directly into the package. Includes a **Tkinter
desktop GUI** for interactive plot and data exploration.

## Quickstart

```bash
pip install -e ".[dev]"           # core + matplotlib + pytest + pyyaml
python examples/basic_usage.py
python -m milne_pn.gui             # launch the desktop GUI
pytest tests/
```

```python
from milne_pn import MilneSolver

solver = MilneSolver(N=41)
print(solver.z0)   # 0.7103982... (Case exact: 0.710446)
```

## The GUI

```bash
python -m milne_pn.gui
# or:  python scripts/run_gui.py
# or, after pip install -e .:  milne-gui
```

A Tkinter desktop app with:

- **Control bar** — pick a PN order and solve instantly; run a convergence
  sweep up to a chosen N_max (runs on a background thread so the window
  stays responsive, with a progress indicator)
- **5 tabs**, each an embedded matplotlib figure with the full pan/zoom/save
  navigation toolbar:
  - *Flux Profile* — φ₀(x), its asymptote, and the extrapolated endpoint
  - *Convergence* — z₀(N) and the log-scale error, after a sweep
  - *Gibbs Ringing* — the boundary angular flux ψ(0,μ) at several orders,
    plus the overshoot-vs-N plateau
  - *Runge Phenomenon* — the classical textbook analogue, side by side
  - *Data Table* — a sortable table of every (N, z₀, error, #modes) row
    from the last convergence sweep
- **File menu** — export all current results as CSVs, save the active plot
  (PNG/PDF/SVG), or quit

Needs `tkinter`, which ships with most Python installs but is a separate OS
package on some minimal Linux distros (`apt install python3-tk` on
Debian/Ubuntu — there's no pip package for it). `tests/unit/test_gui.py`
exercises the GUI logic with a real (or Xvfb virtual) display and skips
itself cleanly if none is available.

## Project layout

```text
milne-pn/
├── src/milne_pn/
│   ├── algorithms/        # PNSystem, MilneSolver, Runge/Gibbs analysis
│   ├── benchmarks/         # convergence + timing studies
│   ├── visualization/      # convergence / flux / Gibbs-ringing plots (matplotlib -> files)
│   ├── gui/                 # Tkinter desktop app (matplotlib -> embedded Tk canvas)
│   └── utils/               # quadrature, config, CSV I/O (numpy-only)
├── tests/{unit,integration,performance}/
├── examples/               # basic_usage.py, benchmark_example.py
├── notebooks/              # algorithm_analysis.ipynb, performance_analysis.ipynb
├── scripts/                # benchmark.py, generate_data.py, run_gui.py
├── data/{input,output}/
├── docs/                   # algorithms.md, architecture.md, performance.md
├── config/config.yaml
├── .github/workflows/tests.yml
├── pyproject.toml, requirements.txt
└── README.md
```

See `docs/architecture.md` for the full rationale behind this layout.

## The physics, briefly

Expanding the angular flux ψ(x,μ) in Legendre moments up to odd order N
turns the transport equation into a linear ODE system `AΦ'(x) + DΦ(x) = 0`.
Its bounded (x→+∞) solutions are an exact analytic asymptotic
(diffusion) mode plus a set of decaying exponential modes, found via
`numpy.linalg.eig`. **Marshak boundary conditions** at x=0 close the
system (`numpy.linalg.solve`), and **z₀ falls out directly as one of the
solved amplitudes** — not from a numerical fit to the computed flux curve.

Full derivation: `docs/algorithms.md`.

## Why convergence is slow: Runge/Gibbs phenomenon

The *exact* boundary angular flux ψ(0,μ) has a jump discontinuity at
μ=0 (zero for μ>0, vacuum; nonzero for μ<0). Representing that
discontinuity with a finite *global* Legendre series is exactly the
setting of **Gibbs' phenomenon** — the reconstruction rings near μ=0 with
an overshoot that does **not** shrink as N grows, only narrows. This
package measures that directly:

```python
from milne_pn.algorithms.milne_solver import MilneSolver
from milne_pn.algorithms.runge_analysis import gibbs_overshoot

for N in (3, 9, 21, 41, 61, 79):
    print(N, gibbs_overshoot(MilneSolver(N))["relative_overshoot"])
# 0.304  0.219  0.203  0.195  0.191  0.191   <- plateau, not decay to 0
```

alongside the classical, purely-mathematical **Runge phenomenon**
(equispaced-node polynomial interpolation of a smooth, peaked function)
as a side-by-side textbook analogue — visible both in
`src/milne_pn/algorithms/runge_analysis.py` /
`notebooks/algorithm_analysis.ipynb` and interactively in the GUI's
*Gibbs Ringing* / *Runge Phenomenon* tabs.

Both are manifestations of the same underlying pathology: a finite-order
**global** polynomial cannot represent a non-smooth function without
oscillating, and that oscillation amplitude is essentially independent of
degree. This is why increasing PN order alone buys only algebraic
(≈O(1/N)), not exponential, convergence of z₀.

## Typical output

```
$ python scripts/generate_data.py
    N         z0 (PN)           error      #modes
    1      0.66666667     -0.04377933         0
    3      0.70509500     -0.00535100         1
    9      0.70964523     -0.00080077         4
   39      0.71039338     -0.00005262        19
   79      0.71043264     -0.00001336        39
```

## Running things

```bash
# core solve, from Python
python examples/basic_usage.py

# convergence + Runge/Gibbs plots
python examples/benchmark_example.py          # -> examples/output/*.png

# interactive desktop GUI
python -m milne_pn.gui

# config-aware CLI tools (read config/config.yaml, override with flags)
python scripts/generate_data.py --nmax 99 --nplot 41
python scripts/benchmark.py --timing

# tests
pytest tests/unit --ignore=tests/unit/test_gui.py
pytest tests/unit/test_gui.py                 # needs a display (or Xvfb)
pytest tests/integration
pytest tests/performance
pytest tests/ --ignore=tests/unit/test_gui.py --cov=milne_pn

# notebooks
jupyter notebook notebooks/algorithm_analysis.ipynb
jupyter notebook notebooks/performance_analysis.ipynb
```

## Reference

The expected converged value, z₀ = 0.7104 mean free paths, is the classical
result derived analytically by K. M. Case via the singular-eigenfunction
expansion of transport theory (Case & Zweifel, *Linear Transport Theory*).
