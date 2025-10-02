"""Polarization abstractions for active microwave scattering."""

from __future__ import annotations

from enum import Enum
from typing import Iterable, Sequence, Tuple

__all__ = [
    "PolarizationState",
    "default_polarization_order",
    "normalize_polarization",
    "normalize_polarization_sequence",
]


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


_DEFAULT_ORDER: Tuple[PolarizationState, ...] = (
    PolarizationState.HH,
    PolarizationState.VV,
    PolarizationState.HV,
    PolarizationState.VH,
)


def default_polarization_order() -> Tuple[PolarizationState, ...]:
    """Return the default ordered tuple of polarization channels."""

    return _DEFAULT_ORDER


def _coerce_polarization(value) -> PolarizationState:
    if isinstance(value, PolarizationState):
        return value
    if isinstance(value, str):
        try:
            return PolarizationState(value.lower())
        except ValueError as exc:  # pragma: no cover - defensive branch
            raise ValueError(f"Unknown polarization string '{value}'") from exc
    raise TypeError(f"Unsupported polarization value: {value!r}")


def normalize_polarization(polarization=None) -> Tuple[PolarizationState, ...]:
    """Normalize a polarization input to an ordered tuple.

    ``polarization`` may be ``None`` (meaning all channels), a
    :class:`PolarizationState`, a string, or an iterable of either.
    """

    if polarization is None:
        return default_polarization_order()

    if isinstance(polarization, (PolarizationState, str)):
        return (_coerce_polarization(polarization),)

    if isinstance(polarization, Iterable):
        return normalize_polarization_sequence(polarization)

    raise TypeError("Polarization must be None, a string, PolarizationState, or an iterable thereof")


def normalize_polarization_sequence(values: Sequence) -> Tuple[PolarizationState, ...]:
    """Normalize an iterable of polarization identifiers preserving order.

    Duplicate entries are removed while preserving the first occurrence.
    """

    seen = set()
    result = []
    for value in values:
        pol = _coerce_polarization(value)
        if pol not in seen:
            seen.add(pol)
            result.append(pol)

    if not result:
        raise ValueError("Polarization sequence must not be empty")

    return tuple(result)
