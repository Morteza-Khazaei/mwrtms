"""Correlation functions used by surface roughness models."""

from __future__ import annotations

from typing import Protocol

import numpy as np

__all__ = [
    "CorrelationFunction",
    "GaussianCorrelation",
    "ExponentialCorrelation",
    "PowerLawCorrelation",
]


class CorrelationFunction(Protocol):
    """Protocol describing roughness correlation functions.

    The strategy pattern is achieved via duck typing: any class implementing
    :meth:`spectrum` can be consumed by :class:`~mwrtms.interface.roughness.SurfaceRoughness`.
    """

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        """Return the *n*-th order roughness spectrum."""
        ...


class _BaseCorrelation:
    """Common validation utilities for concrete correlation functions."""

    __slots__ = ("_ell",)

    def __init__(self, correlation_length: float) -> None:
        if correlation_length <= 0:
            raise ValueError("correlation_length must be positive")
        self._ell = float(correlation_length)

    @property
    def correlation_length(self) -> float:
        return self._ell

    def _validate_order(self, n: int) -> None:
        if n <= 0:
            raise ValueError("order n must be positive")


class GaussianCorrelation(_BaseCorrelation):
    """Gaussian correlation: ρ(τ) = exp(-(τ/ℓ)^2)."""

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        self._validate_order(n)
        K = self._ell * np.hypot(kx, ky)
        return (self._ell**2 / (2.0 * n)) * float(np.exp(-(K**2) / (4.0 * n)))


class ExponentialCorrelation(_BaseCorrelation):
    """Exponential correlation: ρ(τ) = exp(-|τ|/ℓ)."""

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        self._validate_order(n)
        K = self._ell * np.hypot(kx, ky)
        return (self._ell / n) ** 2 * (1.0 + (K / n) ** 2) ** (-1.5)


class PowerLawCorrelation(_BaseCorrelation):
    """Power-law correlation: ρ(τ) = (1 + |τ|/ℓ)^(-p)."""

    __slots__ = ("_p",)

    def __init__(self, correlation_length: float, power: float = 1.5) -> None:
        super().__init__(correlation_length)
        if power <= 0:
            raise ValueError("power must be positive")
        self._p = float(power)

    @property
    def power(self) -> float:
        return self._p

    def spectrum(self, n: int, kx: float, ky: float) -> float:
        self._validate_order(n)
        K = self._ell * np.hypot(kx, ky)
        return (self._ell**2) / (1.0 + (K / max(n, 1)) ** (self._p + 1.0))
