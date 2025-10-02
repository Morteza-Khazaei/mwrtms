"""Exponential correlation function."""

from __future__ import annotations

import numpy as np

from .base import CorrelationFunction

__all__ = ["ExponentialCorrelation"]


class ExponentialCorrelation:
    """Exponential surface correlation ``ρ(τ) = exp(-τ/ℓ)``."""

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        K = np.hypot(kx, ky)
        ell = float(correlation_length)
        n_eff = max(int(n), 1)
        return (ell / n_eff) ** 2 * (1.0 + (ell * K / n_eff) ** 2) ** (-1.5)
