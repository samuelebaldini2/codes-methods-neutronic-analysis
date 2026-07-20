6<!--
SPDX-FileCopyrightText: 2026 Alberto P

SPDX-License-Identifier: CC-BY-NC-SA-4.0
-->

# Instructions

Please consider the following exercises: a corresponsing description is provided.

For each of them create a dedicated `git` branch.
The corresponding development should be composed by multiple commits, in which the contribution of all the members of the working group should be properly attributed.

# Exercises

- [ ] [Software architecture: Mie scattering](#mie-scattering-fortran-wrapping)
- [ ] [Milne problem](#milne-problem-purely-absorbing-medium)
- [ ] [Level-Symmetric quadrature](#level-symmetric-nodes-and-weigths)
- [ ] [2-D diamond finite differences $S_N$](#2-d-diamond-finite-differences-discrete-ordinates)

## Mie scattering: Fortran wrapping
Consider the [MIEV0 wrapper](https://github.com/alberto743/miev0-wrapper), available also on PyPI as package [miescat](https://pypi.org/project/miescat/).
This project is based on the idea of wrapping the original MIEV0 Fortran routine developed by Warren Wiscombe in a single executable that takes the input parameters, such as complex refractive index, grain side, and light wavelength as command line paramters and provide the Legendre moments of the phase function as JSON dictionary in the standard output.
This is in turn wrapped into a convenient Python `subprocess` call, that mimics the behavior of an ordinary Python function.
The objective of this task is to rewrite the wrapping so that the Python processing does not rely on external executables, but is instead based direct function calls, without the need to spawn a new process at each invocation.
This enhanced wrapper may be based on [f2py](https://numpy.org/doc/stable/f2py/) bindings in Fortran or on [pybind11](https://github.com/pybind/pybind11) for C++, via intermediate use of the C interoperability of modern Fortran.
Choose one binding approach; a comparison between f2py and pybind11 is optional/bonus.
`MIEV0` should not contain non-reentrant instructions: it is required to preliminary verify anyhow the presence of `SAVE` or `COMMON` statements, and report back in case before additional actions.
A comparitive analysis of advantages and disadvantages is required, in terms of both effectiveness and performances, e.g. overhead, Mie calculation per unit time with `%timeit`.
The current repository contains some test workflows (GitHub actions) that shall pass. Adaption of relevant definitions are also required, specifically in the CMake definitions.
If it is judged convenient, it is possible to change the build system from `setuptools` to others, e.g. `scikit-build-core` or others.

## Gauss-Legendre quadrature in the unit sphere
Employ the relevant `numpy` and `scipy` routines to construct a product quadrature in the unit sphere, i.e. Gauss-Legendre in the polar angle and Gauss-Chebyshev or trapezoidal in the azimuthal angle.
For a given set of quadrature orders, demonstrate numerically the capability of exact integration of the spherical harmonics, up to a chosen tolerance.
Consider carefully the normalization of the spherical harmonics introduced by the `scipy` library.

## 1-D diamond-difference discrete ordinates
Consider a 1-D slab for which the following initial set of parameters are given:
- length: $d = 5. \, \text{cm}$
- total cross section: $\Sigma_t = 2. \, \text{cm}^{-1}$
- scattering ratio: $c = 0.9$
- external source: $s_\text{ext} = 10 \, \text{cm}^{-3}\text{s}^{-1}$

As a starting parameter, employ a Gauss-Legendre quadrature of order $N = 4$.

Develop a $\text{S}_\text{n}$ solver that is capable of determing the flux in the slab for a source positioned uniformly in the slab or peaked in the center (i.e. uniform in the central region of the mesh chosen, width of your choice to be stated).
Use a diamond-difference spatial scheme with a user-defined number of cells (verify grid independence). Iterate via standard source iteration. As a convergence cryterion, compare the scalar flux in each mesh cell or the right-hand side evaluation of the transport equation, at your choice.

Observe the difference with void conditions in the left and right boundaries, with both reflective conditions, and with mixed void (left) and reflective (right) boundaries.

Modify the scattering of the left and right half so that the scattering ratio is $0.5$ in the left and $0.99$ in the right and repeat the analysis.

Proceed with a linear anisotropic scattering with $g=0.2$ and $p(\mu_0) = \frac{1}{2} (1 + 3 g \mu_0)$. The integral of the scattering cross section shall be equivalent to the $c = 0.9$ case.
Consider also a forward peaked source in the shape of $s_\text{ext}(\mu) = 10 \cdot \frac{1}{2} (1 + \mu) \, \text{cm}^{-3}\text{s}^{-1}$.
Plot the given scattering and source function and describe them.
Update the $\text{S}_\text{n}$ solver to deal with this case as well.

In all cases, observe the effects of the scattering ratio and boundary conditions in the number of iterations.

## Milne problem: purely scattering medium (optional)
Develop a numerical solution of the Milne problem with the $\text{P}_{n}$ (Legendre series expansion) approach for a purely scattering medium.
Use the more convenient `numpy` built-in routines, when applicable.
Employ Marshak-style boundary conditions for the closure of the numerical system.
Observe the convergence of the extrapolation distance (condition at which the scalar flux would go to zero in case of a linear extrapolation, i.e. $0.7104$) increasing the (impair) order of the Legendre polynomial expansion.
