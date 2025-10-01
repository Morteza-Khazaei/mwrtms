"""Base class for scattering mechanisms."""

from __future__ import annotations

from abc import ABC, abstractmethod

from ..core.geometry import ScatteringGeometry
from ..core.polarization import PolarizationState
from ..core.wave import ElectromagneticWave
from ..medium.base import Medium

__all__ = ["ScatteringMechanism"]


class ScatteringMechanism(ABC):
    """Abstract interface for scattering mechanisms implementing a template method."""

    __slots__ = ("_wave", "_geometry", "_last_result")
    _mutable_slots = {"_last_result"}

    def __setattr__(self, name, value):
        if name in self.__slots__ and name not in self._mutable_slots and hasattr(self, name):
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, wave: ElectromagneticWave, geometry: ScatteringGeometry) -> None:
        self._wave = wave
        self._geometry = geometry
        self._last_result: float | None = None

    @property
    def wave(self) -> ElectromagneticWave:
        return self._wave

    @property
    def geometry(self) -> ScatteringGeometry:
        return self._geometry

    def compute(self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState) -> float:
        self._validate_media(medium_above, medium_below)
        result = self._compute_scattering(medium_above, medium_below, polarization)
        self._last_result = max(result, 0.0)
        return self._last_result

    @abstractmethod
    def _compute_scattering(
        self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState
    ) -> float:
        """Subclasses provide the scattering implementation."""

    def _validate_media(self, medium_above: Medium, medium_below: Medium) -> None:
        for medium in (medium_above, medium_below):
            if not isinstance(medium, Medium):
                raise TypeError("media must derive from Medium")
