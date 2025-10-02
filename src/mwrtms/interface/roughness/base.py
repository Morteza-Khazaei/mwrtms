"""Base classes for surface roughness descriptions."""

from __future__ import annotations

from abc import ABC, abstractmethod

__all__ = ["SurfaceRoughness"]


class SurfaceRoughness(ABC):
    """Abstract representation of a surface roughness descriptor."""

    @abstractmethod
    def rms_height(self) -> float:
        """Return the RMS height σ (m)."""

    @abstractmethod
    def correlation_length(self, azimuth_deg: float = 0.0) -> float:
        """Return the correlation length ℓ (m) for the given azimuth."""

    @abstractmethod
    def spectrum(self, n: int, kx: float, ky: float) -> float:
        """Return the roughness spectrum ``W^{(n)}(k_x, k_y)``."""

    @abstractmethod
    def is_isotropic(self) -> bool:
        """Return ``True`` if the roughness is isotropic."""
