"""Surface scattering base classes demonstrating inheritance."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ...core.polarization import PolarizationState
from ...interface.fresnel import FresnelCoefficients
from ...interface.roughness import SurfaceRoughness
from ...medium.base import Medium
from ..base import ScatteringMechanism

__all__ = ["SurfaceScattering"]


class SurfaceScattering(ScatteringMechanism, ABC):
    """Template for surface scattering models."""

    __slots__ = ("_roughness",)
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, wave, geometry, surface_roughness: SurfaceRoughness) -> None:
        super().__init__(wave, geometry)
        self._roughness = surface_roughness

    @property
    def roughness(self) -> SurfaceRoughness:
        return self._roughness

    def _compute_scattering(self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState) -> float:
        fresnel = self._compute_fresnel(medium_above, medium_below)
        kirchhoff = self._compute_kirchhoff(fresnel, polarization)
        complementary = self._compute_complementary(fresnel, polarization)
        return kirchhoff + complementary

    def _compute_fresnel(self, medium_above: Medium, medium_below: Medium) -> dict[str, complex]:
        eps1 = medium_above.permittivity(self.wave.frequency_hz)
        eps2 = medium_below.permittivity(self.wave.frequency_hz)
        theta = self.geometry.theta_i
        return {
            "h": FresnelCoefficients.reflection(eps1, eps2, theta, "h"),
            "v": FresnelCoefficients.reflection(eps1, eps2, theta, "v"),
        }

    @abstractmethod
    def _compute_kirchhoff(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        """Kirchhoff component implemented by subclasses."""

    @abstractmethod
    def _compute_complementary(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        """Complementary component implemented by subclasses."""

    def _surface_spectrum(self, order: int, delta_kx: float, delta_ky: float) -> float:
        return self._roughness.spectrum(order, delta_kx, delta_ky)

    def _prefactor(self) -> float:
        k = self.wave.wavenumber
        theta = self.geometry.theta_i
        return (k**2 * np.cos(theta) ** 4) / (2.0 * np.pi)
