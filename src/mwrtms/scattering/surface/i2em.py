"""Improved IEM implementation."""

from __future__ import annotations

import math

from ...core.polarization import PolarizationState
from ...interface.roughness import SurfaceRoughness
from .base import SurfaceScattering

__all__ = ["I2EMModel"]


class I2EMModel(SurfaceScattering):
    """Improved IEM surface scattering model."""

    MODEL_NAME = "I2EM"
    __slots__ = ("_transition_alpha",)
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, wave, geometry, surface_roughness: SurfaceRoughness, transition_alpha: float = 3.0) -> None:
        super().__init__(wave, geometry, surface_roughness)
        self._transition_alpha = float(transition_alpha)

    def _compute_kirchhoff(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        base = super()._prefactor()
        ks = self.wave.wavenumber * self.roughness.rms_height
        transition = 1.0 / (1.0 + math.exp(-self._transition_alpha * (ks - 1.2)))
        pol_key = "v" if polarization in (PolarizationState.VV, PolarizationState.V) else "h"
        return base * abs(fresnel[pol_key]) ** 2 * transition

    def _compute_complementary(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        if polarization in (PolarizationState.HV, PolarizationState.VH):
            return 0.0
        ks = self.wave.wavenumber * self.roughness.rms_height
        spectrum = self._surface_spectrum(1, 0.0, 0.0)
        return 0.1 * ks**2 * spectrum * abs(fresnel["v"] - fresnel["h"]) ** 2
