"""Power-law correlation function."""

from __future__ import annotations

import numpy as np
from scipy.special import kv, gamma

from .base import CorrelationFunction

__all__ = ["PowerLawCorrelation"]


class PowerLawCorrelation:
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
