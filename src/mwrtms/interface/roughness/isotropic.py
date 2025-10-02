"""Isotropic surface roughness."""

from __future__ import annotations

from dataclasses import dataclass

from ..correlation.base import CorrelationFunction
from .base import SurfaceRoughness

__all__ = ["IsotropicRoughness"]


@dataclass(frozen=True)
class IsotropicRoughness(SurfaceRoughness):
    """Encapsulate isotropic surface roughness parameters."""

    rms_height_m: float
    correlation_length_m: float
    correlation_function: CorrelationFunction

    def __post_init__(self) -> None:
        if self.rms_height_m <= 0.0 or self.correlation_length_m <= 0.0:
            raise ValueError("rms_height_m and correlation_length_m must be positive")

    def rms_height(self) -> float:
        return self.rms_height_m

    def correlation_length(self, azimuth_deg: float = 0.0) -> float:
        return self.correlation_length_m

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        return self.correlation_function.spectrum(n, kx, ky, self.correlation_length_m)

    def is_isotropic(self) -> bool:
        return True
