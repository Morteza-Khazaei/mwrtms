"""Physical constants used throughout the mwRTMs foundation layer."""

from __future__ import annotations

import math

__all__ = ["SPEED_OF_LIGHT", "EPSILON_0", "MU_0", "PLANCK_CONSTANT"]

SPEED_OF_LIGHT: float = 2.998e8
"""Speed of light in vacuum (m/s)."""

EPSILON_0: float = 8.8541878128e-12
"""Vacuum permittivity ε₀ (F/m)."""

MU_0: float = 4.0 * math.pi * 1e-7
"""Vacuum permeability μ₀ (H/m)."""

PLANCK_CONSTANT: float = 6.62607015e-34
"""Planck constant h (J·s)."""
