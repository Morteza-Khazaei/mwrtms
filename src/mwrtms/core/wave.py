"""Fundamental electromagnetic wave abstractions."""

from __future__ import annotations

from typing import Dict, Optional, TYPE_CHECKING

import numpy as np

from .constants import SPEED_OF_LIGHT

if TYPE_CHECKING:  # pragma: no cover
    from .geometry import ScatteringGeometry

__all__ = ["ElectromagneticWave"]


class ElectromagneticWave:
    """Representation of a monochromatic electromagnetic wave.

    The class demonstrates **encapsulation** through private attributes and
    exposes derived properties via read-only accessors. Expensive calculations
    are cached lazily so downstream models can reuse the object efficiently.
    """

    __slots__ = ("_frequency_hz", "_wavelength_cache", "_angular_frequency_cache")
    _mutable_slots = {"_wavelength_cache", "_angular_frequency_cache"}

    def __init__(self, frequency_hz: float) -> None:
        if frequency_hz <= 0:
            raise ValueError("frequency_hz must be positive")
        self._frequency_hz: float = float(frequency_hz)
        self._wavelength_cache: Optional[float] = None
        self._angular_frequency_cache: Optional[float] = None

    def __setattr__(self, name, value):
        if name in self.__slots__ and name not in self._mutable_slots and hasattr(self, name):
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    @property
    def frequency_hz(self) -> float:
        """Carrier frequency in Hertz."""
        return self._frequency_hz

    @property
    def frequency_ghz(self) -> float:
        """Carrier frequency expressed in GHz."""
        return self._frequency_hz / 1e9

    @property
    def angular_frequency(self) -> float:
        """Angular frequency (rad/s)."""
        if self._angular_frequency_cache is None:
            self._angular_frequency_cache = 2.0 * np.pi * self._frequency_hz
        return self._angular_frequency_cache

    @property
    def wavelength(self) -> float:
        """Free-space wavelength (m)."""
        if self._wavelength_cache is None:
            self._wavelength_cache = SPEED_OF_LIGHT / self._frequency_hz
        return self._wavelength_cache

    @property
    def wavenumber(self) -> float:
        """Free-space wavenumber (rad/m)."""
        return 2.0 * np.pi / self.wavelength

    @property
    def band(self) -> str:
        """IEEE radar band designation."""
        f_ghz = self.frequency_ghz
        if 1 <= f_ghz < 2:
            return "L"
        if 2 <= f_ghz < 4:
            return "S"
        if 4 <= f_ghz < 8:
            return "C"
        if 8 <= f_ghz < 12:
            return "X"
        if 12 <= f_ghz < 18:
            return "Ku"
        if 18 <= f_ghz < 27:
            return "K"
        if 27 <= f_ghz < 40:
            return "Ka"
        return "Unknown"

    def wave_vectors(self, geometry: "ScatteringGeometry") -> Dict[str, np.ndarray]:
        """Compute incident and scattered wave vectors."""
        return geometry.wave_vectors(self.wavenumber)
