"""Kirchhoff's law validator."""

from __future__ import annotations

from ..result.brightness import BrightnessTemperatureResult
from ..result.scattering import ScatteringResult

__all__ = ["KirchhoffLawValidator"]


class KirchhoffLawValidator:
    """Validate ε + |R|^2 ≈ 1 for flat surfaces."""

    @staticmethod
    def check_flat_surface(
        scattering_result: ScatteringResult,
        emissivity_result: BrightnessTemperatureResult,
        tolerance: float = 0.01,
    ) -> bool:
        hv = ["hh", "vv"]
        for pol in hv:
            if pol not in scattering_result.data:
                continue
            sigma0 = scattering_result[pol]
            theta = scattering_result.geometry.theta_i
            reflection_estimate = sigma0 / (4.0 * 3.141592653589793 * (abs(theta) + 1e-12))
            if pol == "hh":
                emissivity = emissivity_result.emissivity_h
            else:
                emissivity = emissivity_result.emissivity_v
            if abs((emissivity + reflection_estimate) - 1.0) > tolerance:
                return False
        return True
