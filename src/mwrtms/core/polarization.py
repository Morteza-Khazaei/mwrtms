"""Polarization abstractions."""

from __future__ import annotations

from enum import Enum
from typing import Dict

__all__ = ["PolarizationState", "StokesVector"]


class PolarizationState(Enum):
    """Standard active and passive microwave polarizations."""

    HH = "hh"
    VV = "vv"
    HV = "hv"
    VH = "vh"
    H = "h"
    V = "v"

    @property
    def is_copol(self) -> bool:
        return self in (PolarizationState.HH, PolarizationState.VV)

    @property
    def is_crosspol(self) -> bool:
        return self in (PolarizationState.HV, PolarizationState.VH)

    @property
    def is_passive(self) -> bool:
        return self in (PolarizationState.H, PolarizationState.V)


class StokesVector:
    """Stokes vector representation with encapsulated state."""

    __slots__ = ("_i", "_q", "_u", "_v")

    def __init__(self, i: float, q: float, u: float, v: float) -> None:
        self._i = float(i)
        self._q = float(q)
        self._u = float(u)
        self._v = float(v)

    @property
    def I(self) -> float:  # noqa: N802 - Physical notation
        return self._i

    @property
    def Q(self) -> float:  # noqa: N802 - Physical notation
        return self._q

    @property
    def U(self) -> float:  # noqa: N802 - Physical notation
        return self._u

    @property
    def V(self) -> float:  # noqa: N802 - Physical notation
        return self._v

    @property
    def brightness_h(self) -> float:
        return (self._i + self._q) / 2.0

    @property
    def brightness_v(self) -> float:
        return (self._i - self._q) / 2.0

    def to_tb(self) -> Dict[str, float]:
        return {"TBh": self.brightness_h, "TBv": self.brightness_v}
