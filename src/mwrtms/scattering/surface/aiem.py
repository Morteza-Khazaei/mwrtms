"""Advanced Integral Equation Method (AIEM) surface model."""

from __future__ import annotations

import math
import numpy as np

from ...factory import register_model
from .base import SurfaceScattering

__all__ = ["AIEMModel"]


@register_model("aiem")
class AIEMModel(SurfaceScattering):
    """AIEM implementation with simplified complementary term."""

    MODEL_NAME = "AIEM"

    def __init__(self, wave, geometry, surface_roughness, *, n_max: int = 15) -> None:
        super().__init__(wave, geometry, surface_roughness)
        if n_max <= 0:
            raise ValueError("n_max must be positive")
        self._n_max = int(n_max)

    def _compute_kirchhoff(self, R_h, R_v, polarization):
        theta = self._geometry.theta_i_rad
        k = self._wave.wavenumber
        sigma = self._roughness.rms_height()

        if polarization.value in ("vv", "v"):
            R = R_v
        elif polarization.value in ("hh", "h"):
            R = R_h
        else:
            return 0.0

        prefactor = (k**2 * np.cos(theta) ** 4) / (2.0 * np.pi)
        slope_prob = np.exp(-(k * sigma * np.cos(theta)) ** 2)

        kx = -2.0 * k * np.sin(theta)
        ky = 0.0
        spectrum = self._roughness.spectrum(0, kx, ky)

        return prefactor * abs(R) ** 2 * slope_prob * spectrum

    def _compute_complementary(self, R_h, R_v, polarization):
        k = self._wave.wavenumber
        theta = self._geometry.theta_i_rad
        sigma = self._roughness.rms_height()

        total = 0.0
        for n in range(1, min(self._n_max, 10) + 1):
            kx = -2.0 * k * np.sin(theta)
            ky = 0.0
            W_n = self._roughness.spectrum(n, kx, ky)
            I_n = 0.1 * (sigma * k) ** n
            factorial_n = math.factorial(n)
            total += (sigma**2 / factorial_n) * abs(I_n) ** 2 * W_n

        return (k**2 / (8.0 * np.pi)) * total
