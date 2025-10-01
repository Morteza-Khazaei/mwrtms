"""Base medium abstraction."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

__all__ = ["Medium"]


class Medium(ABC):
    """Abstract electromagnetic medium."""

    @abstractmethod
    def permittivity(self, frequency_hz: float, temperature_k: float | None = None) -> complex:
        raise NotImplementedError

    @abstractmethod
    def permeability(self, frequency_hz: float) -> complex:
        raise NotImplementedError

    @abstractmethod
    def physical_temperature(self) -> float:
        raise NotImplementedError

    def impedance(self, frequency_hz: float) -> complex:
        eps = self.permittivity(frequency_hz)
        mu = self.permeability(frequency_hz)
        return np.sqrt(mu / eps)

    def loss_tangent(self, frequency_hz: float, temperature_k: float | None = None) -> float:
        eps = self.permittivity(frequency_hz, temperature_k)
        return abs(eps.imag / eps.real)
