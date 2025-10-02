"""Abstract scattering mechanism base classes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

__all__ = ["ScatteringMechanism"]

if TYPE_CHECKING:
    from ..core import PolarizationState
    from ..core.radar_modes import RadarConfiguration
    from ..medium import Medium


class ScatteringMechanism(ABC):
    """Universal base class for scattering mechanisms."""

    def __init__(self, wave, geometry) -> None:
        self._wave = wave
        self._geometry = geometry

    @abstractmethod
    def compute(self, medium_above, medium_below, polarization) -> float:
        """Return the backscatter coefficient σ⁰ (linear power)."""

        raise NotImplementedError

    def compute_with_config(
        self,
        medium_above: "Medium",
        medium_below: "Medium",
        polarization: "PolarizationState",
        config: "RadarConfiguration",
    ) -> float:
        """Compute scattering using an explicit radar configuration."""

        if not config.matches_geometry(self._geometry):
            raise ValueError(
                "Radar configuration geometry does not match the scattering mechanism geometry"
            )
        return self.compute(medium_above, medium_below, polarization)

    @property
    def wave(self):
        return self._wave

    @property
    def geometry(self):
        return self._geometry
