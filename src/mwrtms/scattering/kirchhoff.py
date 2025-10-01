"""Kirchhoff approximation implementation."""

from __future__ import annotations

import math

import numpy as np

from ..core.geometry import ScatteringGeometry
from ..core.polarization import PolarizationState
from ..core.wave import ElectromagneticWave
from ..interface.fresnel import FresnelCoefficients
from ..interface.roughness import SurfaceRoughness
from ..medium.base import Medium
from .base import ScatteringMechanism

__all__ = ["KirchhoffApproximation"]


class KirchhoffApproximation(ScatteringMechanism):
    """Kirchhoff (geometric optics) surface scattering."""

    model_name = "kirchhoff"

    def compute(
        self,
        wave: ElectromagneticWave,
        geometry: ScatteringGeometry,
        interface: SurfaceRoughness,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
    ) -> float:
        k = wave.wavenumber
        theta = geometry.theta_i
        eps1 = medium_above.permittivity(wave.frequency_hz)
        eps2 = medium_below.permittivity(wave.frequency_hz)
        pol = "h" if polarization in (PolarizationState.HH, PolarizationState.HV, PolarizationState.H) else "v"
        R = FresnelCoefficients.reflection(eps1, eps2, theta, pol)
        sigma = interface.rms_height
        cos_theta = math.cos(theta)
        prefactor = (k**2 * cos_theta**4) / (2.0 * math.pi)
        slope_probability = math.exp(-((k * sigma * cos_theta) ** 2))
        spec = interface.spectrum(1, -2.0 * k * math.sin(theta), 0.0)
        return float(prefactor * abs(R) ** 2 * slope_probability * spec)
