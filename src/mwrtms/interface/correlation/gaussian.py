"""Gaussian correlation function."""

from __future__ import annotations

import numpy as np

from .base import CorrelationFunction

__all__ = ["GaussianCorrelation"]


class GaussianCorrelation:
    """Gaussian surface correlation ``ρ(τ) = exp(-τ²/ℓ²)``."""

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        K = np.hypot(kx, ky)
        ell = float(correlation_length)
        n_eff = max(int(n), 1)
        return (ell ** 2 / (2.0 * n_eff)) * np.exp(-ell ** 2 * K ** 2 / (4.0 * n_eff))
