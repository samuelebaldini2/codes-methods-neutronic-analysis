# Algorithms

## The Milne problem

A semi-infinite, purely-scattering half-space (x in [0, infinity), scattering
ratio c = 1, no absorption, no internal source), with the source
effectively at x -> infinity. In mean-free-path units, the 1D transport equation
(azimuthally symmetric, isotropic scattering) is:

```
mu d(psi)/dx + psi(x,mu) = (1/2) phi_0(x)
```

where phi_0(x) = integral_{-1}^{1} psi(x,mu) dmu is the scalar flux. We want
the bounded (as x->infinity, up to a linear asymptote) solution subject to a
vacuum boundary at x=0: no incoming radiation, psi(0,mu) = 0 for mu > 0.

## PN expansion

Expand the angular flux in Legendre moments:

```
psi(x,mu) = sum_{l=0}^{N} (2l+1)/2 * phi_l(x) * P_l(mu)        (N odd)
```

Projecting the transport equation onto each P_l(mu) and using the recursion
mu*P_l = [(l+1)P_{l+1} + l*P_{l-1}]/(2l+1) gives the PN moment equations:

```
l = 0:      phi_1(x) = 0
l = 1..N:   l/(2l+1)*phi_{l-1}(x) + (l+1)/(2l+1)*phi_{l+1}(x) + phi_l(x) = 0
```

with closure phi_{N+1} = 0. Written as a first-order linear system:

```
A Phi(x) + D Phi(x) = 0
```

A is tridiagonal with zero diagonal (A[l,l-1] = l/(2l+1), A[l,l+1] =
(l+1)/(2l+1)), and D = diag(0,1,1,...,1).
(algorithms/pn_system.py :: PNSystem)

## Exponential modes and the eigenproblem

Substituting Phi(x) = v * exp(-omega*x) turns the ODE system into the
generalized eigenproblem D v = omega A v, i.e. an ordinary eigenproblem for
M = A^{-1}D (np.linalg.solve + np.linalg.eig,
algorithms/milne_solver.py :: MilneSolver).

The double eigenvalue omega = 0 is defective (a genuine 2x2 Jordan block,
not diagonalizable) -- no numerical eigensolver reliably extracts it. Its
exact closed form is used directly instead:

```
phi_0(x) = P0 + q0*x,   phi_1(x) = -q0/3,   phi_l(x) = 0  for l >= 2
```

-- the diffusion asymptote, exact within any PN truncation. All other
(nonzero) eigenvalues are real for this problem (verified numerically:
max_imag_residual is 0 at every tested order) and occur in +/- pairs; only
those with Re(omega) > 0 are physical, since the solution must stay bounded
as x -> +infinity.

## Marshak boundary conditions

Enforcing psi(0,mu) = 0 pointwise for all mu > 0 is impossible with a finite
Legendre truncation. Marshak closure instead enforces it weakly, in the
odd-order half-range moments:

```
integral_0^1 P_k(mu) psi(0,mu) dmu = 0,   k = 1, 3, 5, ..., N
```

using the half-range overlap matrix H[k,l] = integral_0^1 P_k(mu)P_l(mu)dmu,
computed exactly via Gauss-Legendre quadrature (numpy.polynomial.legendre.leggauss,
utils/quadrature.py). This gives exactly (N+1)/2 linear equations for the
(N+1)/2 unknowns: P0 and the amplitudes of the (N-1)/2 decaying discrete
modes -- solved with np.linalg.solve.

## z0 is read off directly

With q0 normalized to 1, the asymptotic flux is phi_0(x) = x + P0, whose
zero-crossing is at x = -P0. So z0 = P0, one of the amplitudes returned
by the linear solve -- an analytic construction, not a numerical fit to
computed flux values.

## Convergence, and why it is slow: the Runge/Gibbs phenomenon

algorithms/runge_analysis.py addresses this directly. The exact
boundary angular flux has a jump discontinuity at mu = 0: identically
zero for mu > 0, finite and nonzero for mu < 0. The PN method reconstructs
psi(0,mu) as a finite-degree global Legendre series -- approximating a
discontinuous function with a global polynomial basis is precisely the
setting of Gibbs phenomenon: the reconstruction rings near the
discontinuity, and -- crucially -- the ring amplitude does not shrink
with N, only its width does.

This is the angular-domain sibling of the classical Runge phenomenon
(polynomial interpolation at equispaced nodes ringing near the domain
edges as degree grows): both are pathologies of representing a non-smooth
or steep function with a single global polynomial, and both explain why
increasing PN order alone buys only algebraic (approx O(1/N)), not
exponential, convergence.

gibbs_overshoot() quantifies this: the relative overshoot in the (should
be exactly zero) mu > 0 region measured at N = 3, 9, 21, 41, 61, 79 is
30%, 22%, 20%, 19.5%, 19.1%, 19.1% -- a plateau, not a decay to
zero. runge_interpolation_demo() reproduces the classical textbook case
(interpolating 1/(1+25x^2) at equispaced nodes) side by side, for direct
visual comparison -- see visualization/angular_flux_visualizer.py,
notebooks/algorithm_analysis.ipynb, and the GUI Gibbs Ringing / Runge
Phenomenon tabs.

## Reference

The expected converged value, z0 = 0.7104 mean free paths, is the classical
result derived analytically by K. M. Case via the singular-eigenfunction
expansion of transport theory (Case & Zweifel, Linear Transport Theory).
