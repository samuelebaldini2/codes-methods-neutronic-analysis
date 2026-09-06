# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

import os
import subprocess
import json
from pathlib import Path
from typing import Tuple, Union
import numpy as np


exec_name = "miescat.exe" if os.name == "nt" else "miescat"
DEFAULT_EXECUTABLE = Path(__file__).parent / exec_name


def run_console(program: Path, *args: str) -> Tuple[str, str, int]:
    result = subprocess.run(
        [str(program)] + list(args),
        capture_output=True,
        text=True,
        check=True
    )
    return result.stdout, result.stderr, result.returncode


def compute_mie_scattering(m_real: float,
                           m_img: float,
                           radius: float,
                           wavelength: float,
                           executable: Union[Path, str]=DEFAULT_EXECUTABLE
                           ) -> dict:
    executable_path = Path(executable)
    output, _, _ = run_console(
        executable_path,
        str(m_real),
        str(m_img),
        str(radius),
        str(wavelength),
    )
    return json.loads(output)


class Mie:
    '''Mie scattering calculator.'''

    def __init__(self, m_real: float, m_img: float, radius: float, wavelength: float):
        self._m_real = m_real
        self._m_img = m_img
        self._radius = radius
        self._wavelength = wavelength
        self._results = compute_mie_scattering(m_real, m_img, radius, wavelength)
        self._legendre_moments = None

    @property
    def refractive_index(self):
        return complex(self._m_real, -self._m_img)

    @property
    def size_parameter(self):
        return 2 * np.pi * self._radius / self._wavelength

    @property
    def radius(self):
        return self._radius

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def sphere_geometric_cross_section(self):
        return np.pi * self._radius ** 2

    @property
    def extinction_efficiency(self):
        return self._results['extinction_efficiency']

    @property
    def scattering_efficiency(self):
        return self._results['scattering_efficiency']

    @property
    def absorption_efficiency(self):
        return self._results['absorption_efficiency']

    @property
    def asymmetry_factor(self):
        return self._results['asymmetry_factor']

    @property
    def radiation_pressure_efficiency(self):
        return self._results['radiation_pressure_efficiency']

    @property
    def single_scattering_albedo(self):
        return self._results['single_scattering_albedo']

    @property
    def legendre_moments(self):
        if self._legendre_moments is None:
            self._legendre_moments = np.array(self._results['legendre_moments'])
        return self._legendre_moments
