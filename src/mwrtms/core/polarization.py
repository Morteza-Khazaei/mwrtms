"""Polarization abstractions for active microwave scattering."""

from __future__ import annotations

from enum import Enum

__all__ = ["PolarizationState"]


class PolarizationState(Enum):
    """Enumerate co- and cross-polarised radar channels."""

    HH = "hh"
    VV = "vv"
    HV = "hv"
    VH = "vh"

    @property
    def is_copol(self) -> bool:
        """Return ``True`` for co-polarised channels (HH, VV)."""

        return self in (PolarizationState.HH, PolarizationState.VV)

    @property
    def is_crosspol(self) -> bool:
        """Return ``True`` for cross-polarised channels (HV, VH)."""

        return self in (PolarizationState.HV, PolarizationState.VH)
