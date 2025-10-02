"""Dobson et al. (1985) soil dielectric mixing model."""

from __future__ import annotations

import numpy as np

from .base import DielectricModel

__all__ = ["DobsonModel"]


class DobsonModel(DielectricModel):
    """Dobson et al. (1985) dielectric mixing model for moist soils."""

    def compute(
        self,
        frequency_hz: float,
        *,
        moisture: float,
        clay: float,
        sand: float,
        temperature_k: float = 293.15,
    ) -> complex:
        """Return the complex permittivity for the supplied soil state."""

        rho_b = 1.3 + 0.6 * clay  # bulk density (g/cm³)
        rho_s = 2.65  # particle density (g/cm³)

        eps_winf = 4.9
        eps_w0 = 88.045 - 0.4147 * (temperature_k - 273.15)
        tau_w = (
            1.1109e-10
            - 3.824e-12 * (temperature_k - 273.15)
            + 6.938e-14 * (temperature_k - 273.15) ** 2
        )

        omega = 2 * np.pi * frequency_hz
        eps_fw = eps_winf + (eps_w0 - eps_winf) / (1 + 1j * omega * tau_w)

        alpha = 0.65
        eps_s = 4.7 + 0j

        eps_soil = (
            1
            + (rho_b / rho_s) * (eps_s**alpha - 1)
            + moisture * eps_fw**alpha
            - moisture
        ) ** (1 / alpha)

        return eps_soil
