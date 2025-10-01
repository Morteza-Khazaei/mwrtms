"""Frequency and wavenumber helper classes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import SPEED_OF_LIGHT

__all__ = ["FrequencyBand", "Wavenumber"]


@dataclass(frozen=True)
class FrequencyBand:
    """Spectral band definition."""

    start_hz: float
    stop_hz: float

    def __post_init__(self) -> None:
        if self.start_hz <= 0 or self.stop_hz <= 0:
            raise ValueError("Frequencies must be positive")
        if self.stop_hz <= self.start_hz:
            raise ValueError("stop_hz must be greater than start_hz")

    @property
    def center_hz(self) -> float:
        return 0.5 * (self.start_hz + self.stop_hz)

    @property
    def bandwidth(self) -> float:
        return self.stop_hz - self.start_hz

    def contains(self, frequency_hz: float) -> bool:
        return self.start_hz <= frequency_hz <= self.stop_hz


@dataclass(frozen=True)
class Wavenumber:
    """Helper converting between frequency and wavenumber."""

    value: float

    @classmethod
    def from_frequency(cls, frequency_hz: float) -> "Wavenumber":
        return cls(2.0 * np.pi * frequency_hz / SPEED_OF_LIGHT)

    @property
    def wavelength(self) -> float:
        return 2.0 * np.pi / self.value
