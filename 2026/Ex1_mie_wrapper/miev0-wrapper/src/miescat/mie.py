# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

import os
import subprocess
import json
from pathlib import Path
from typing import Tuple, Union
import numpy as np

from . import _mie_wrapper


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
                           executable: Union[Path, str] = DEFAULT_EXECUTABLE
                           ) -> dict:

    output = _mie_wrapper.compute_mie(
        m_real,
        m_img,
        radius,
        wavelength,
    )

    number_of_moments = int(output[8])
    legendre_moments = output[9][:number_of_moments + 1].tolist()

    results = {
        "refractive_index_real": float(np.float32(m_real)),
        "refractive_index_imaginary": float(np.float32(m_img)),
        "particle_radius": float(np.float32(radius)),
        "wavelength": float(np.float32(wavelength)),
        "size_parameter": float(output[0]),
        "extinction_efficiency": float(output[1]),
        "scattering_efficiency": float(output[2]),
        "absorption_efficiency": float(output[3]),
        "sphere_geometric_cross_section": float(output[4]),
        "single_scattering_albedo": float(output[5]),
        "asymmetry_factor": float(output[6]),
        "radiation_pressure_efficiency": float(output[7]),
        "number_of_moments": number_of_moments,
        "legendre_moments": legendre_moments,
    }

    return results


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
