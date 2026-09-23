"""
milne_pn
--------
A P_N (spherical-harmonics / Legendre) solver for the Milne problem: a
semi-infinite, purely-scattering half-space with the source at x -> +inf.
"""

__version__ = "0.1.0"

from .algorithms.milne_solver import MilneSolver
from .algorithms.pn_system import PNSystem

__all__ = ["MilneSolver", "PNSystem", "__version__"]
