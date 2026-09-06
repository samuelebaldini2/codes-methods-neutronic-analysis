# SPDX-FileCopyrightText: 2026 Alberto P
#
# SPDX-License-Identifier: MPL-2.0

import math
import numpy as np
import pytest
import miescat


def test_default_executable_exists():
    assert miescat.default_executable.exists(), f"Executable not found at {miescat.default_executable}"
    assert miescat.default_executable.is_file()
    assert miescat.default_executable.name in {"miescat", "miescat.exe"}


def test_mie_class_computes_scattering():
    mieclass = miescat.Mie(1.33, 0.01, 1.0, 0.55)

    assert isinstance(mieclass.refractive_index, complex)
    assert mieclass.radius == pytest.approx(1.0)
    assert mieclass.wavelength == pytest.approx(0.55)
    assert mieclass.size_parameter == pytest.approx(2 * math.pi * 1.0 / 0.55)
    assert mieclass.extinction_efficiency >= 0
    assert mieclass.scattering_efficiency >= 0
    assert mieclass.absorption_efficiency >= 0
    assert 0 <= mieclass.single_scattering_albedo <= 1
    assert isinstance(mieclass.legendre_moments, np.ndarray)
    assert mieclass.legendre_moments.ndim == 1
    assert mieclass.legendre_moments.size > 0


def test_compute_mie_scattering_function():
    mieres = miescat.compute_mie_scattering(1.3484, 0.001, 1.0, 0.41)

    assert isinstance(mieres, dict)
    assert len(mieres['legendre_moments']) == 55
    assert 0 <= mieres['single_scattering_albedo'] <= 1
