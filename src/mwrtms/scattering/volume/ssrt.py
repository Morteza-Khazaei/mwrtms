"""Single-scattering radiative transfer (SSRT) model."""

from __future__ import annotations

import math

from ...core.polarization import PolarizationState
from ...medium.base import Medium
from .base import VolumeScattering

__all__ = ["SSRTModel"]


class SSRTModel(VolumeScattering):
    """Simplified SSRT implementation coupling soil and vegetation."""

    MODEL_NAME = "SSRT"
    __slots__ = ("_canopy_tau", "_single_scatter_albedo", "_soil_backscatter")
    _mutable_slots: set[str] = set()

    def __init__(self, wave, geometry, canopy_optical_depth: float, single_scatter_albedo: float, soil_backscatter: dict[str, float]) -> None:
        super().__init__(wave, geometry)
        if canopy_optical_depth < 0:
            raise ValueError("canopy_optical_depth must be non-negative")
        if not 0 <= single_scatter_albedo <= 1:
            raise ValueError("single_scatter_albedo must be in [0, 1]")
        self._canopy_tau = float(canopy_optical_depth)
        self._single_scatter_albedo = float(single_scatter_albedo)
        self._soil_backscatter = {k.lower(): float(v) for k, v in soil_backscatter.items()}

    def _compute_volume_response(self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState) -> float:
        pol_key = polarization.value
        soil_term = self._soil_backscatter.get(pol_key, 0.0)
        attenuation = math.exp(-2.0 * self._canopy_tau / max(math.cos(self.geometry.theta_i), 1e-3))
        volume = self._single_scatter_albedo * (1.0 - attenuation)
        return soil_term * attenuation + volume
