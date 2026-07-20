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

## Milne problem: purely absorbing medium


## Level-Symmetric nodes and weigths


## 2-D diamond finite differences discrete ordinates
