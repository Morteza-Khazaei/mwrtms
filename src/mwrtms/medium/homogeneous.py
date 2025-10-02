"""Simple homogeneous medium helpers."""

from __future__ import annotations

from .base import Medium

__all__ = ["HomogeneousMedium"]


class HomogeneousMedium(Medium):
    """Medium with constant complex permittivity."""

    __slots__ = ("_permittivity",)

    def __init__(self, permittivity: complex, temperature_k: float = 293.15) -> None:
        super().__init__(temperature_k)
        self._permittivity = complex(permittivity)

    def permittivity(self, frequency_hz: float) -> complex:  # pragma: no cover - trivial accessor
        return self._permittivity
