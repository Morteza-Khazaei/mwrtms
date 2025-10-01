"""Flat surface emissivity calculation."""

from __future__ import annotations

from ..interface.fresnel import FresnelCoefficients

__all__ = ["FlatSurfaceEmissivity"]


class FlatSurfaceEmissivity:
    """Compute emissivity of a planar interface using Fresnel coefficients."""

    @staticmethod
    def compute(eps1: complex, eps2: complex, theta: float, polarization: str) -> float:
        return FresnelCoefficients.emissivity(eps1, eps2, theta, polarization)
