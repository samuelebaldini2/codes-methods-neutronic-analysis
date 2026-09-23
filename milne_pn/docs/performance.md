# Performance

## What dominates the cost

`MilneSolver(N)` does three dense-linear-algebra operations on
(N+1)x(N+1) matrices:

1. `np.linalg.solve(A, D)` to form `M = A^{-1}D`
2. `np.linalg.eig(M)` for the full eigendecomposition (the expensive
   step: a general, non-symmetric dense eigenproblem)
3. `np.linalg.solve(Mtx, rhs)` for the Marshak closure -- much smaller,
   only (N+1)/2 unknowns

Plus a Gauss-Legendre quadrature of size `4x(N+1)` for the half-range
overlap matrix `H` (`PNSystem`), which is cheap relative to the above.

## Measured timings

(Single run, sandbox CPU; `repeats=3`, mean +/- std. Reproduce with
`python scripts/benchmark.py --timing` or
`milne_pn.benchmarks.benchmark.timing_benchmark`.)

| N   | mean time | std   |
|-----|-----------|-------|
| 11  | 2.4 ms    | 1.7 ms|
| 31  | 5.3 ms    | 0.3 ms|
| 61  | 17.0 ms   | 0.2 ms|
| 101 | 44.8 ms   | 0.2 ms|
| 151 | 113.6 ms  | 5.2 ms|

Fitting `time ~ c*N^p` over this range gives an empirical exponent
p ~ 1.5. The asymptotic complexity of dense `solve`/`eig` on an nxn
matrix is O(n^3); the measured exponent is well below that here because
these matrix sizes (<=152x152) are still small enough that fixed overhead
and LAPACK own practical (better-than-worst-case) behavior dominate over
the asymptotic cubic term. Do not extrapolate p~1.5 to much larger N.

In absolute terms: the full convergence study used throughout this
project (N = 1, 3, ..., 79, forty separate solves) completes in well
under a second, and even N=151 is ~0.1 s. Performance is not a practical
constraint for this project purpose -- including inside the GUI, where a
single-N solve (the `Solve` button) is effectively instantaneous and only
the multi-order convergence sweep is run on a background thread (purely
so the window stays responsive/redraws during it, not because it is
actually slow).

## Where NOT to over-optimize

  * `PNSystem.H` (the Marshak overlap matrix) is recomputed independently
    for each N via quadrature. Cheap relative to the eigendecomposition.
  * The Runge/Gibbs analysis (`runge_analysis.py`) reconstructs
    `psi(0,mu)` on a fine `n_mu`-point grid (default 2000, 1000 in the
    GUI Gibbs tab for snappier redraws) for each PN order inspected.
    Cheap O(n_mu * N) work, not a bottleneck.

## If you need to go faster

The natural next step for much larger N (hundreds+) would be to exploit
the tridiagonal structure of `A` and `D` directly (a proper generalized
eigenvalue solver for banded/tridiagonal pencils, e.g. via `scipy.linalg`
banded routines) instead of forming the dense `M = A^{-1}D` and calling
the general-matrix `np.linalg.eig`. This project intentionally keeps to
dense numpy built-ins for simplicity and clarity, since the PN orders
that actually matter for this problem (up to ~100) are already
comfortably fast without it.
