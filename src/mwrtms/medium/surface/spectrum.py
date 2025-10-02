"""Spectral helpers for :class:`~mwrtms.medium.surface.base.Surface`."""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .base import Surface


def density(surface: Surface) -> np.ndarray:
    """Return the power spectral density (PSD)."""

    return surface.power_spectrum()


def radial(surface: Surface) -> Tuple[np.ndarray, np.ndarray]:
    """Return azimuthally averaged PSD."""

    return surface.radial_spectrum()


def frequency_grid(surface: Surface) -> Tuple[np.ndarray, np.ndarray]:
    """Return the angular frequency grid for the surface."""

    return surface.frequency_grid()
