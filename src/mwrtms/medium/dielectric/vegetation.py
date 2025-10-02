"""Vegetation dielectric material model."""

from __future__ import annotations

import numpy as np

from .base import DielectricModel

__all__ = ["VegetationMaterialModel"]


class VegetationMaterialModel(DielectricModel):
    """Empirical vegetation dielectric model."""

    def compute(
        self,
        frequency_hz: float,
        *,
        gravimetric_moisture: float,
        temperature_k: float = 293.15,
    ) -> complex:
        """Return the complex permittivity for vegetation material."""

        f_ghz = frequency_hz / 1e9
        mg = gravimetric_moisture

        eps_real = (1.7 + 3.2 * mg + 6.5 * mg**2) * (1 - 0.002 * (temperature_k - 293.15))
        eps_imag = (0.7 * mg * f_ghz / (1 + (f_ghz / 18) ** 2)) * (1 + 0.003 * (temperature_k - 293.15))

        return complex(eps_real, eps_imag)
