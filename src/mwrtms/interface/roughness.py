"""Surface roughness abstractions."""

from __future__ import annotations

from typing import Dict

from .correlation import CorrelationFunction

__all__ = ["SurfaceRoughness"]


class SurfaceRoughness:
    """Surface roughness parameters with encapsulated state."""

    __slots__ = ("_rms_height", "_correlation_length", "_correlation")
    _mutable_slots: set[str] = set()

    def __init__(
        self,
        rms_height: float,
        correlation_length: float,
        correlation_function: CorrelationFunction,
    ) -> None:
        if rms_height < 0:
            raise ValueError("rms_height must be non-negative")
        if correlation_length <= 0:
            raise ValueError("correlation_length must be positive")
        self._rms_height = float(rms_height)
        self._correlation_length = float(correlation_length)
        self._correlation = correlation_function

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    @property
    def rms_height(self) -> float:
        return self._rms_height

    @property
    def correlation_length(self) -> float:
        return self._correlation_length

    @property
    def correlation_function(self) -> CorrelationFunction:
        return self._correlation

    def normalized_parameters(self, wavenumber: float) -> Dict[str, float]:
        return {
            "ks": wavenumber * self._rms_height,
            "kl": wavenumber * self._correlation_length,
        }

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        return self._correlation.spectrum(n, kx, ky)

    def is_smooth(self, wavenumber: float, threshold: float = 0.3) -> bool:
        return wavenumber * self._rms_height < threshold

    def is_rough(self, wavenumber: float, threshold: float = 3.0) -> bool:
        return wavenumber * self._rms_height > threshold
