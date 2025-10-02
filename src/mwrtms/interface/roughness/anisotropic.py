"""Anisotropic surface roughness."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..correlation.base import CorrelationFunction
from .base import SurfaceRoughness

__all__ = ["AnisotropicRoughness"]


@dataclass(frozen=True)
class AnisotropicRoughness(SurfaceRoughness):
    """Directional roughness suitable for tilled soils."""

    rms_height_m: float
    correlation_length_x_m: float
    correlation_length_y_m: float
    correlation_function: CorrelationFunction

    def __post_init__(self) -> None:
        if self.rms_height_m <= 0.0:
            raise ValueError("rms_height_m must be positive")
        if self.correlation_length_x_m <= 0.0 or self.correlation_length_y_m <= 0.0:
            raise ValueError("correlation lengths must be positive")

    def rms_height(self) -> float:
        return self.rms_height_m

    def correlation_length(self, azimuth_deg: float = 0.0) -> float:
        phi = np.deg2rad(azimuth_deg)
        cos2 = np.cos(phi) ** 2
        sin2 = np.sin(phi) ** 2
        return 1.0 / (cos2 / self.correlation_length_x_m + sin2 / self.correlation_length_y_m)

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        ell = self.correlation_length(np.degrees(np.arctan2(ky, kx)))
        return self.correlation_function.spectrum(n, kx, ky, ell)

    def is_isotropic(self) -> bool:
        return False

    @property
    def anisotropy_ratio(self) -> float:
        return max(self.correlation_length_x_m, self.correlation_length_y_m) / min(
            self.correlation_length_x_m, self.correlation_length_y_m
        )
