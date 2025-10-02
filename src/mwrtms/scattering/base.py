"""Abstract scattering mechanism base classes."""

from __future__ import annotations

from abc import ABC, abstractmethod

__all__ = ["ScatteringMechanism"]


class ScatteringMechanism(ABC):
    """Universal base class for scattering mechanisms."""

    def __init__(self, wave, geometry) -> None:
        self._wave = wave
        self._geometry = geometry

    @abstractmethod
    def compute(self, medium_above, medium_below, polarization) -> float:
        """Return the backscatter coefficient σ⁰ (linear power)."""

        raise NotImplementedError

    @property
    def wave(self):
        return self._wave

    @property
    def geometry(self):
        return self._geometry
