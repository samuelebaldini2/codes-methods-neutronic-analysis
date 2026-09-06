<!--
SPDX-FileCopyrightText: 2026 Alberto P

SPDX-License-Identifier: MPL-2.0
-->

# miescat - MIEV0 wrapper

Fortran and Python wrapper of the MIEV0 code to compute Mie scattering by Warren J. Wiscombe.

The original MIEV0 code has been retrieved from `ftp://climate1.gsfc.nasa.gov/wiscombe/Single_Scatt/Homogen_Sphere/Exact_Mie/`.

## Fortran executable
The Fortran executable is named `miescat`.
It computes the Mie scattering parameters and the Legendre moments starting from a real and imaginary part of the scattering frequency, a radius, and a wavelength.
The output may be parsed as json with [jq](https://github.com/jqlang/jq).

A Python simple wrapper script `mieleg` may be used to have more explicit command line argugments.
```bash
mieleg --m_real 1.3484 --m_img 0.001 --radius 1.0 --wavelength 0.41
```

## Python package
The Python package is a wrapper around the Fortran executable and exposes a Python interface to compute Mie scattering.
The execution is spawn as a subprocess that the resuting Mie coefficients and Legendre moments are returned as a Python dictionary.

The package installs a command-line entry point `miescat`.
It exposes the `compute_mie_scattering` function and the `Mie` class.

The package is available on PyPI [miescat](https://pypi.org/project/miescat/).
```shell
pip install miescat
```

## Requirements
- Python 3.9+
- CMake
- Fortran compiler

## License
- The present wrapping tool is released under the terms of the Mozilla Public License 2.0 (MPL 2.0) license.
- MIEV0 is distributed by NASA and apparently belongs to the public domain. All credit goes to Warren J. Wiscombe.
