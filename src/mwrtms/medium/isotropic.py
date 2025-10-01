"""Concrete isotropic medium implementation."""

from __future__ import annotations

from typing import Optional

from .base import Medium

__all__ = ["IsotropicMedium"]


class IsotropicMedium(Medium):
    """Simple isotropic medium with encapsulated properties."""

    __slots__ = ("_eps", "_mu", "_temperature")
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, permittivity: complex, temperature_k: float, permeability: complex = 1.0) -> None:
        self._eps = complex(permittivity)
        self._mu = complex(permeability)
        if temperature_k <= 0:
            raise ValueError("temperature_k must be positive")
        self._temperature = float(temperature_k)

    def permittivity(self, frequency_hz: float, temperature_k: Optional[float] = None) -> complex:
        return self._eps

    def permeability(self, frequency_hz: float) -> complex:
        return self._mu

    def physical_temperature(self) -> float:
        return self._temperature
