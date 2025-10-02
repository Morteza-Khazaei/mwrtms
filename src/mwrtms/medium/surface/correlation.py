"""Autocorrelation function helpers for surface statistics."""

from __future__ import annotations

from typing import Protocol

import numpy as np
from scipy.special import kv, gamma

__all__ = [
    "CorrelationFunction",
    "Exponential",
    "Gaussian",
    "PowerLaw",
]


class CorrelationFunction(Protocol):
    """Protocol for surface correlation functions."""

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        """Return the roughness spectrum ``W^{(n)}(k_x, k_y)``."""

        ...


class Exponential:
    """Exponential surface correlation ``ρ(τ) = exp(-τ/ℓ)``."""

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        K = np.hypot(kx, ky)
        ell = float(correlation_length)
        n_eff = max(int(n), 1)
        return (ell / n_eff) ** 2 * (1.0 + (ell * K / n_eff) ** 2) ** (-1.5)


class Gaussian:
    """Gaussian surface correlation ``ρ(τ) = exp(-τ²/ℓ²)``."""

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        K = np.hypot(kx, ky)
        ell = float(correlation_length)
        n_eff = max(int(n), 1)
        return (ell ** 2 / (2.0 * n_eff)) * np.exp(-ell ** 2 * K ** 2 / (4.0 * n_eff))


class PowerLaw:
    """Power-law correlation ``ρ(τ) = (1 + |τ|/ℓ)^{-p}``."""

    def __init__(self, power: float = 1.5) -> None:
        if power <= 0:
            raise ValueError("power must be positive")
        self._power = float(power)

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        K = np.hypot(kx, ky)
        if K == 0.0:
            return float("inf")
        n_eff = max(int(n), 1)
        exponent = self._power * n_eff - 1.0
        try:
            bessel = (correlation_length * K / 2.0) ** exponent * kv(-exponent, correlation_length * K)
            return (correlation_length ** 2 / gamma(self._power * n_eff)) * bessel
        except Exception:
            return 0.0
