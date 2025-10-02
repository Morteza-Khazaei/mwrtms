"""Base class for surface scattering models."""

from __future__ import annotations

from abc import abstractmethod

from ..base import ScatteringMechanism

__all__ = ["SurfaceScattering"]


class SurfaceScattering(ScatteringMechanism):
    """Template method for surface scattering models."""

    def __init__(self, wave, geometry, surface_roughness) -> None:
        super().__init__(wave, geometry)
        self._roughness = surface_roughness

    def compute(self, medium_above, medium_below, polarization) -> float:
        R_h, R_v = self._compute_fresnel(medium_above, medium_below)
        kirchhoff = self._compute_kirchhoff(R_h, R_v, polarization)
        complementary = self._compute_complementary(R_h, R_v, polarization)
        return kirchhoff + complementary

    @abstractmethod
    def _compute_kirchhoff(self, R_h, R_v, polarization):
        """Return the Kirchhoff contribution."""

    @abstractmethod
    def _compute_complementary(self, R_h, R_v, polarization):
        """Return the complementary (multiple) contribution."""

    def _compute_fresnel(self, medium_above, medium_below):
        from ...interface import FresnelCoefficients

        eps1 = medium_above.permittivity(self._wave.frequency_hz)
        eps2 = medium_below.permittivity(self._wave.frequency_hz)
        theta = self._geometry.theta_i_rad

        R_h = FresnelCoefficients.reflection(eps1, eps2, theta, "h")
        R_v = FresnelCoefficients.reflection(eps1, eps2, theta, "v")
        return R_h, R_v
