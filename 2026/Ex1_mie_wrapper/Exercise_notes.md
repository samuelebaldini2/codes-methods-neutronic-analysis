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
## Validation

After implementing the new wrapper, the existing test suite was executed with:

```bash
pytest test/
```

and all tests passed successfully. This verifies that the new implementation preserves the expected public interface and output structure.
The numerical comparison was then performed for the reference case specified in the exercise:

-complex refractive index: m = 1.02239 - 0.0119i

-particle radius: r = 2 µm

-wavelength: λ = 0.5 µm

The same input data were used for both implementations.

## Numerical comparison

The two implementations give very similar numerical results.

For scalar quantities, the absolute and relative differences are small.

| Quantity | Absolute error | Relative error |
|---|---:|---:|
| Size parameter | 2.81e-08 | 1.12e-09 |
| Extinction efficiency | 2.83e-10 | 2.72e-10 |
| Scattering efficiency | 2.18e-10 | 4.46e-10 |
| Absorption efficiency | 4.99e-10 | 9.01e-10 |
| Single-scattering albedo | 2.58e-10 | 5.51e-10 |
| Asymmetry factor | 3.80e-10 | 3.82e-10 |
| Radiation-pressure efficiency | 7.72e-11 | 1.39e-10 |


For the Legendre moments, the comparison gives:

- maximum absolute error: 2.82e-08;
- mean absolute error: 4.59e-09;
- maximum relative error: 5.11e-08;
- mean relative error: 2.01e-08.

These implementations can be considered numerically consistent.
The complete numerical comparison is stored in `numerical_comparison.csv`.

## Performance comparison

To compare the performance of the two approaches, the total execution time was measured for an increasing number of successive calls:

`N = 1, 5, 10, 20`

For each value of `N`, the same Mie calculation was repeated using the old subprocess-based wrapper and the new direct `f2py` binding.

The measured results are:

| Number of calls | Old wrapper total time [s] | New wrapper total time [s] | Old time/call [s] | New time/call [s] | Speedup |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.2325 | 0.2185 | 0.2325 | 0.2185 | 1.06 |
| 5 | 1.0998 | 1.0333 | 0.2200 | 0.2067 | 1.06 |
| 10 | 2.4424 | 2.0326 | 0.2442 | 0.2033 | 1.20 |
| 20 | 4.4533 | 4.0559 | 0.2227 | 0.2028 | 1.10 |

The corresponding calculation rates are approximately:

| Number of calls | Old wrapper [calls/s] | New wrapper [calls/s] |
|---:|---:|---:|
| 1 | 4.30 | 4.58 |
| 5 | 4.55 | 4.84 |
| 10 | 4.09 | 4.92 |
| 20 | 4.49 | 4.93 |

The new wrapper is consistently faster for all tested values of N.

The direct wrapper reduces the per-call overhead because it avoids:

- creation of a new operating-system process;
- conversion of the input parameters to command-line strings;
- communication through standard output.

However, the measured speedup is moderate. In this case, the actual `MIEV0` computation represents a significant fraction of the total execution time, so the subprocess overhead is only part of the cost of each call. As a consequence, removing that overhead produces a measurable but relatively limited reduction of the total execution time.

The full timing results are stored in `timing_comparison.csv`.
