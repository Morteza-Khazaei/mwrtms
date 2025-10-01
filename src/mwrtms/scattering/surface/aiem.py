"""Advanced Integral Equation Method (AIEM) implementation."""

from __future__ import annotations

import math

from ...core.polarization import PolarizationState
from ...interface.roughness import SurfaceRoughness
from .base import SurfaceScattering

__all__ = ["AIEMModel"]


class AIEMModel(SurfaceScattering):
    """AIEM surface scattering model demonstrating inheritance and encapsulation."""

    MODEL_NAME = "AIEM"
    __slots__ = ("_include_multiple", "_series_order")
    _mutable_slots: set[str] = set()

    def __init__(self, wave, geometry, surface_roughness: SurfaceRoughness, include_multiple: bool = True) -> None:
        super().__init__(wave, geometry, surface_roughness)
        self._include_multiple = bool(include_multiple)
        self._series_order = 6

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def _compute_kirchhoff(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        pol_key = "v" if polarization in (PolarizationState.VV, PolarizationState.V) else "h"
        reflection = fresnel[pol_key]
        theta = self.geometry.theta_i
        slope = math.exp(-(self.wave.wavenumber * self.roughness.rms_height * math.cos(theta)) ** 2)
        return self._prefactor() * abs(reflection) ** 2 * slope

    def _compute_complementary(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        if polarization in (PolarizationState.HV, PolarizationState.VH):
            return 0.0
        weight = 0.0
        ks = self.wave.wavenumber * self.roughness.rms_height
        for order in range(1, self._series_order + 1):
            spectrum = self._surface_spectrum(order, 0.0, 0.0)
            weight += (ks ** (2 * order) / math.factorial(order)) * spectrum
        if not self._include_multiple:
            weight *= 0.5
        return 0.25 * abs(fresnel["v"] - fresnel["h"]) ** 2 * weight
