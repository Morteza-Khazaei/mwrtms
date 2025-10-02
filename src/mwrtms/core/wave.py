"""Fundamental electromagnetic wave abstractions."""

from __future__ import annotations

from dataclasses import dataclass
from math import pi

from .constants import SPEED_OF_LIGHT

__all__ = ["ElectromagneticWave"]


@dataclass(frozen=True)
class ElectromagneticWave:
    """Immutable description of a monochromatic electromagnetic wave.

    Encapsulation is achieved by keeping the underlying frequency private and
    exposing read-only derived quantities via properties. Immutable dataclasses
    naturally prevent accidental mutation after construction, which aligns with
    the physical interpretation of a fixed-frequency plane wave.
    """

    frequency_hz: float

    def __post_init__(self) -> None:
        if self.frequency_hz <= 0.0:
            raise ValueError("frequency_hz must be positive")

    @property
    def wavelength(self) -> float:
        """Free-space wavelength λ = c / f [m]."""

        return SPEED_OF_LIGHT / self.frequency_hz

    @property
    def wavenumber(self) -> float:
        """Free-space wavenumber k = 2π / λ [rad/m]."""

        return 2.0 * pi / self.wavelength

    @property
    def frequency_ghz(self) -> float:
        """Frequency expressed in GHz."""

        return self.frequency_hz / 1e9
