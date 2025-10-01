"""Physical constants used throughout the mwRTMs package."""

from __future__ import annotations

import math

__all__ = [
    "SPEED_OF_LIGHT", "BOLTZMANN", "VACUUM_PERMITTIVITY",
    "VACUUM_PERMEABILITY", "PLANCK", "FERMI_TEMP"
]

SPEED_OF_LIGHT: float = 299_792_458.0
"""Speed of light in vacuum (m/s)."""

BOLTZMANN: float = 1.380_649e-23
"""Boltzmann constant (J/K)."""

VACUUM_PERMITTIVITY: float = 8.854_187_8128e-12
"""Permittivity of free space (F/m)."""

VACUUM_PERMEABILITY: float = 4 * math.pi * 1e-7
"""Permeability of free space (H/m)."""

PLANCK: float = 6.626_070_15e-34
"""Planck constant (J·s)."""

FERMI_TEMP: float = 2.725
"""Cosmic microwave background temperature (K)."""
