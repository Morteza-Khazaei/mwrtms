"""Mironov et al. (2009) soil dielectric model."""

from __future__ import annotations

import numpy as np

from .base import DielectricModel

__all__ = ["MironovModel"]


class MironovModel(DielectricModel):
    """Generalised refractive mixing dielectric model for moist soils."""

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

        f_ghz = frequency_hz / 1e9

        eps_winf = 4.9
        eps_w0 = 88.045 - 0.4147 * (temperature_k - 273.15)
        tau_w = 1.1109e-10 - 3.824e-12 * (temperature_k - 273.15)

        omega = 2 * np.pi * frequency_hz
        eps_fw = eps_winf + (eps_w0 - eps_winf) / (1 + 1j * omega * tau_w)

        n_d = 1.634 - 0.539 * clay + 0.2748 * clay**2
        n_soil = n_d + (np.sqrt(eps_fw) - 1) * moisture

        return n_soil**2
