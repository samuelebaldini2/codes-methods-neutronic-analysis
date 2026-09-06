# Exercise 1 - Mie scattering software interface wrapping
## Context of the exercise
We started from a Python interface for a Fortran-based code. This interface created a new processes at runtime, passing the input data as command-line arguments and then reading the standard output, converting it to a python dictionary.
## Aim of the exercise
The aim of the exercise was to modify the wrapping structure by creating a Python binding for the Fortran code, thus avoiding the need to create a new process at each call.
## New wrapper design
The original Fortran program was based on the `MIEV0` routine. The input arguments required by `MIEV0` were derived from the parameters read from the command line, and this preprocessing was performed in the main program.

In order to create a Python binding while preserving the same computational structure, we introduced a new `compute_mie()` routine in `mie_wrapper.F90`. This routine performs the same preliminary operations previously carried out by the main program and then calls `MIEV0`.

A Python binding was then generated for `compute_mie()` using `f2py` and linked to the `src` package. Finally, the `compute_mie_scattering()` function in `mie.py` was modified to call the bound `compute_mie()` routine directly, instead of launching the old executable through a subprocess. The python dictionary is created within `compute_mie_scattering()` from the values returned by `compute_mie()`.

## Project structure
<small>
<pre>
├── Exercise_notes.md
├── main.py
└── miev0-wrapper
    ├── miev
    │   ├── _mie_wrapper.cpython-311-x86_64-linux-gnu.so  ← <strong>compiled Python extension module</strong>
    │   └── miewrapper.F90                                ← <strong>implementation of compute_mie()</strong>
    └── src
        └── miescat
            ├── _mie_wrapper.cpython-311-x86_64-linux-gnu.so
            │   └── <strong>symbolic link</strong> ─────► ../../miev/_mie_wrapper.cpython-311-x86_64-linux-gnu.so
            ├── mie_old.py  ← <strong>old implementation</strong>
            └── mie.py      ← <strong>new implementation using compute_mie()</strong>
</pre>
</small>
## Comparative analysis
