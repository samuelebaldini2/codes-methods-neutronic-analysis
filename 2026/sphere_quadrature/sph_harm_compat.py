"""
Compatibility shim for scipy.special.sph_harm_y.

scipy.special.sph_harm_y was added in SciPy 1.15 and became the only
option once the older scipy.special.sph_harm was removed in SciPy
1.17. Many installations still run older SciPy (this is common with
Anaconda's bundled version) and only have the legacy sph_harm -- which
has a DIFFERENT argument order AND a swapped theta/phi convention:

    sph_harm_y(n, m, theta, phi)   theta = polar,     phi = azimuthal   (SciPy >= 1.15)
    sph_harm(m, n, theta, phi)     theta = azimuthal, phi = polar      (SciPy < 1.17, deprecated)

(confirmed directly against the SciPy docs for both functions).

This module exposes a single function, sph_harm_y(n, m, theta, phi),
always using the modern (n, m, theta=polar, phi=azimuthal) convention
that the rest of this project uses (see PHYSICS_QUADRATURE.md, section
5.1) -- re-mapped onto whichever of the two SciPy functions is
actually available.
"""

try:
    from scipy.special import sph_harm_y

except ImportError:
    from scipy.special import sph_harm as _legacy_sph_harm

    def sph_harm_y(n, m, theta, phi):
        return _legacy_sph_harm(m, n, phi, theta)
