"""Base emission model definitions."""

from __future__ import annotations

from abc import ABC, abstractmethod

from ..core.geometry import ScatteringGeometry
from ..core.polarization import PolarizationState
from ..core.wave import ElectromagneticWave
from ..medium.base import Medium

__all__ = ["EmissionModel"]


class EmissionModel(ABC):
    """Abstract base class for emissivity models."""

    model_name: str = "emission-model"

    @abstractmethod
    def compute_emissivity(
        self,
        wave: ElectromagneticWave,
        geometry: ScatteringGeometry,
        medium: Medium,
        polarization: PolarizationState,
    ) -> float:
        raise NotImplementedError

    def compute_brightness_temperature(self, emissivity: float, physical_temperature: float) -> float:
        return emissivity * physical_temperature
