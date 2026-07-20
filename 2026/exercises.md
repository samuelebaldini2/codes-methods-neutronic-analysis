<!--
SPDX-FileCopyrightText: 2026 Alberto P

SPDX-License-Identifier: CC-BY-NC-SA-4.0
-->

# Instructions

Please consider the following exercises: a corresponding description is provided for each of them.

For each exercise create a dedicated `git` branch.
The corresponding development shall be composed of multiple commits, in which the contribution of all the members of the working group is properly attributed.

**Deliverable.** The deliverable consists exclusively of the `git` branches in the group repository: source code, tests, and any figures or short notes required by the single exercise.
No separate report is expected; explanatory content shall be provided as Markdown files or Jupyter notebooks committed in the branch.

**Deadline.** Third week of September 2026 (soft deadline; an extension may be granted upon request).

**Organisation.**
The four exercises are largely independent and can be developed in parallel on separate branches.
Note that the Mie exercise involves the build toolchain: it is recommended to raise questions if environment-related issues emerge.
You may consider the quadrature exercise and the isotropic part of the discrete ordinates solver as an intermediate milestone.
The working group is composed of 5 members: you may choose your preferred organization.

# Exercises

- [ ] [Mie scattering software interface wrapping](#mie-scattering-software-interface-wrapping)
- [ ] [Gauss-Legendre quadrature in the unit sphere](#gauss-legendre-quadrature-in-the-unit-sphere)
- [ ] [1-D diamond-difference discrete ordinates solver](#1-d-diamond-difference-discrete-ordinates-solver)
- [ ] [Milne problem and extrapolation distance (optional)](#milne-problem-and-extrapolation-distance-optional)


## Mie scattering software interface wrapping

The goal of this exercise is to understand the peculiarities of the different methodologies to manage the interfaces, in the context of the software architecture choices of numerical simulation codes.
Please refer also to the [SATIF16](https://github.com/alberto743/satif16) project for more information.

Consider the [MIEV0 wrapper](https://github.com/alberto743/miev0-wrapper), available also on PyPI as package [miescat](https://pypi.org/project/miescat/).
This project is based on the idea of wrapping the original Mie computation Fortran routine developed by Warren Wiscombe `MIEV0.f` via a single executable `mieleg.F90`, that takes the input parameters, i.e. complex refractive index, grain size, and light wavelength, as command line parameters and provides the Legendre moments of the phase function as a JSON dictionary in the standard output.
This is in turn wrapped into a convenient `subprocess` call in the Python `run_console` function, which is used by `compute_mie_scattering` to mimic the behaviour of an ordinary Python function.
Also, the `Mie` class can be used in the context of object-oriented programming, with a Pythonic style.

The objective of this task is to rewrite the wrapping so that the Python processing does not rely on external executables, but is instead based on direct function calls, without the need to spawn a new process at each invocation.

As a preliminary point to be verified, the underlying function shall be re-entrant: it is required to check the presence of `SAVE` or `COMMON` statements in the `MIEV0` routine.
Should any of them be present, open a GitHub issue before proceeding.
Indeed, the modified wrapper shall not introduce any new internal state.

From the methodological point of view, this enhanced wrapper may be based on [f2py](https://numpy.org/doc/stable/f2py/) bindings in Fortran, or on [pybind11](https://github.com/pybind/pybind11) for C++, via intermediate use of the C interoperability of modern Fortran (i.e. the `iso_c_binding` intrinsic module).
Choose one binding approach; a comparison between `f2py` and `pybind11` is optional/bonus.

The public Python API shall remain unchanged: `compute_mie_scattering` and the `Mie` class shall keep the same signature and return the same data structures, so that the existing test suite remains applicable without modification.

A comparative analysis of advantages and disadvantages between the current approach and the new wrapper is required, in terms of both effectiveness and performance, e.g. overhead, Mie calculations per unit time with `%timeit`.
For this comparison, consider the calculation of the Mie Legendre coefficients for a particle having refractive index $m = 1.02239 - 0.0119j$, radius $2 \, \mu\text{m}$, and wavelength $\lambda = 0.5 \, \mu\text{m}$.
In order to separate the per-call overhead from the actual computational cost, study the scaling of the total execution time as a function of the number of successive calls.

The current repository contains dedicated test workflows (GitHub Actions, e.g. `.github/workflows/ci.yml`), based on the execution of appropriate `CTest` and `pytest` instructions for testing the wrapper.
The workflows shall be able to run the same test suite with the new wrapper.
Adaptation of the relevant definitions may be required, specifically in the CMake definitions, given the different wrapping strategy.
If it is judged more convenient, it is possible to change the build system from `setuptools` to others, e.g. `scikit-build-core`.


## Gauss-Legendre quadrature in the unit sphere

Employ the relevant `numpy` and `scipy` routines to construct a product quadrature in the unit sphere, i.e. Gauss-Legendre in the polar angle and Gauss-Chebyshev or trapezoidal in the azimuthal angle.

For a given set of quadrature orders, i.e. $N_\mu$ and $N_\phi$, demonstrate numerically the capability of exact integration of the spherical harmonics, up to a chosen tolerance.
For instance, you shall be able to demonstrate the possibility to integrate exactly over $\mu$ the polynomial orders up to $(2 N_\mu - 1)$.
Consider carefully the normalization of the spherical harmonics introduced by the `scipy` library, i.e. by looking at the documentation of the `sph_harm_y` function.

As a bonus, demonstrate mathematically the equivalence between Gauss-Chebyshev and trapezoidal integration in the azimuthal angle.


## 1-D diamond-difference discrete ordinates solver

Consider the transport equation for a 1-D slab geometry.

The following initial description of the system is given:

- length: $d = 5. \, \text{cm}$
- total cross section: $\Sigma_t = 2. \, \text{cm}^{-1}$
- scattering ratio: $c = 0.9$
- external source: $s_\text{ext} = 10 \, \text{cm}^{-3}\text{s}^{-1}$

The goal of this exercise is to construct an $\text{S}_\text{N}$ solver based on a diamond-difference spatial scheme.
As a starting setting, employ a Gauss-Legendre quadrature of order $N = 4$ for the angular variable $\mu$.

Two sets of sources shall be considered:

1. uniformly distributed sources in the whole slab;
2. peaked sources at the centre of the slab.

Since the solver is based on a finite-difference scheme, the delta-like source is to be represented as uniform over a central region of the mesh, with a width of your choice, to be stated.
The total integral source strength shall be preserved with respect to the uniform case, so that the two configurations are directly comparable.

Three sets of boundary conditions are foreseen:

1. void conditions on the left and right boundaries, i.e. void/void;
2. reflective conditions on both sides, i.e. reflect/reflect;
3. mixed void (left) and reflective (right) boundaries, i.e. void/reflect.

Develop an $\text{S}_\text{N}$ solver that is capable of determining the flux in the slab, using the most convenient `numpy` features in Python.
Use a diamond-difference spatial scheme with a user-defined number of cells (verify grid independence).
Iterate via standard source iteration.

As a convergence criterion, compare the angular flux moments $\varphi_l$ in each mesh cell between successive iterations, or alternatively the right-hand side evaluation of the transport equation.
Use a relative tolerance of $10^{-8}$ in the infinity norm as a default value, so that iteration counts are comparable across configurations.
A comparison between the two methodologies is bonus/optional.

The diamond-difference scheme is not positivity-preserving: negative cell-edge fluxes may appear when the optical thickness of a cell is large, in particular with the peaked source.
This is an expected feature of the scheme, not a coding error.
In such case, reduce the cell size until the effect disappears; the appearance of negative values is itself a useful indicator for mesh refinement.

In order to prove the correctness of the solver, verify that the solution for the uniform source with reflect/reflect boundary conditions is equal to the analytical reference solution $\varphi_0 = \frac{s_\text{ext}}{\Sigma_t (1 - c)}$, where $\varphi_0$ is the scalar flux.

The user shall be able to represent graphically the solution for the scalar flux.
Observe the shape of the scalar flux when changing the source definition and the boundary conditions.

Modify the scattering of the left and right half so that the scattering ratio is $0.5$ on the left and $0.99$ on the right, and repeat the analysis, preserving the total cross section.

By modifying $c$ from $0.2$ to $0.99$, plot the number of iterations required as a function of the scattering ratio, for the three sets of boundary conditions, and comment on the differences.

### Anisotropic case

Proceed by considering a linear anisotropic extension of the solver that has been developed, namely:

- linear anisotropic scattering with $g = 0.2$ and $p(\mu_0) = \frac{1}{2} (1 + 3 g \mu_0)$;
- forward peaked source in the shape of $s_\text{ext}(\mu) = 10 \cdot \frac{1}{2} (1 + \mu) \, \text{cm}^{-3}\text{s}^{-1}$.

The scattering kernel is expanded as

$$\Sigma_s(\mu_0) = \sum_{l} \frac{2l+1}{2} \, \Sigma_{s,l} \, P_l(\mu_0),$$

with $\Sigma_{s,0} = c \, \Sigma_t = 1.8 \, \text{cm}^{-1}$ and $\Sigma_{s,1} = g \, \Sigma_{s,0} = 0.36 \, \text{cm}^{-1}$, so that the integral of the scattering cross section is equivalent to the $c = 0.9$ case.
The same $\frac{2l+1}{2}$ factor shall appear in the reconstruction of the scattering source from the flux moments $\varphi_l$.

Note that $\int_{-1}^{1} s_\text{ext}(\mu) \, d\mu = 10 \, \text{cm}^{-3}\text{s}^{-1}$: the zeroth moment is preserved with respect to the isotropic case, so that the two configurations are directly comparable.
Only the first moment is introduced.

Plot the given scattering and source functions and describe them.
Update the $\text{S}_\text{N}$ solver to deal with this case as well and observe the results graphically.


## Milne problem and extrapolation distance (optional)

Develop a numerical solution of the Milne problem (semi-infinite half space with the source at $x \to +\infty$) with the $\text{P}_N$ (Legendre series expansion) approach, for a purely scattering medium.

The $\text{P}_N$ equations in the half space admit an analytical solution in $x$ as a superposition of exponential modes $e^{-\omega_i x}$, where the $\omega_i$ are obtained from the eigenvalues of the $\text{P}_N$ coefficient matrix; the physical solution retains only the modes that do not diverge as $x \to +\infty$, plus the linear asymptotic term.
Employ Marshak-style boundary conditions at $x = 0$ for the closure of the numerical system and for fixing the mode amplitudes.
Use the most convenient `numpy` built-in routines, when applicable (e.g. for the eigenvalue problem and for the solution of the resulting linear system).

The extrapolation distance $z_0$ is the intercept with $\varphi_0 = 0$ of the linear extrapolation of the asymptotic part of the scalar flux, and is to be obtained directly from the analytical expression rather than from a numerical fit.

Observe the convergence of $z_0$ when increasing the (odd) order of the Legendre polynomial expansion.
The expected value is $0.7104$ mean free paths, as demonstrated analytically by Case.
Note that convergence is slow: several $\text{P}_N$ orders may be needed.
