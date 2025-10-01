"""Classical Integral Equation Method (IEM)."""

from __future__ import annotations

from ...core.polarization import PolarizationState
from ...interface.roughness import SurfaceRoughness
from .base import SurfaceScattering

__all__ = ["IEMModel"]


class IEMModel(SurfaceScattering):
    """Classical IEM model (baseline implementation)."""

    MODEL_NAME = "IEM"
    __slots__ = ()

    def _compute_kirchhoff(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        pol_key = "v" if polarization in (PolarizationState.VV, PolarizationState.V) else "h"
        return self._prefactor() * abs(fresnel[pol_key]) ** 2

    def _compute_complementary(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        spectrum = self._surface_spectrum(1, 0.0, 0.0)
        if polarization in (PolarizationState.HV, PolarizationState.VH):
            return 0.05 * spectrum
        return 0.2 * spectrum * abs(fresnel["v"] - fresnel["h"]) ** 2
