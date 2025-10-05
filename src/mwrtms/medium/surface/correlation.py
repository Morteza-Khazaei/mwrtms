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
    """Exponential surface correlation ``ρ(τ) = exp(-τ/ℓ)``.
    
    For anisotropic surfaces, the correlation length varies with azimuth:
    L(φ) = Lx cos²(φ) + Ly sin²(φ)
    
    where φ is the azimuthal angle in the spectral domain.
    
    References
    ----------
    Yang, Y., & Chen, K. S. (2019). Polarized backscattering from spatially
    anisotropic rough surface. IEEE TGRS, 57(9), 6608-6618.
    Equation (1) and (4).
    """

    def __init__(
        self,
        correlation_length_x: float | None = None,
        correlation_length_y: float | None = None,
    ) -> None:
        """Initialize exponential correlation function.
        
        Parameters
        ----------
        correlation_length_x : float, optional
            Correlation length along x-axis (minor axis). If None, uses isotropic.
        correlation_length_y : float, optional
            Correlation length along y-axis (major axis). If None, uses isotropic.
        """
        self._lx = correlation_length_x
        self._ly = correlation_length_y
        self._is_anisotropic = (
            correlation_length_x is not None
            and correlation_length_y is not None
            and not np.isclose(correlation_length_x, correlation_length_y)
        )

    def spectrum(
        self,
        n: int,
        kx: float,
        ky: float,
        correlation_length: float,
    ) -> float:
        """Compute nth-order roughness spectrum W^(n)(kx, ky).
        
        Parameters
        ----------
        n : int
            Spectral order (1, 2, 3, ...)
        kx : float
            Spatial frequency x-component
        ky : float
            Spatial frequency y-component
        correlation_length : float
            Isotropic correlation length (used if not anisotropic)
            
        Returns
        -------
        float
            Roughness spectrum value
            
        Notes
        -----
        For isotropic surfaces:
            W^(n)(K) = (L/n)² [1 + (KL/n)²]^(-1.5)
            
        For anisotropic surfaces:
            W^(n)(K,φ) = (L(φ)/n)² [1 + (KL(φ)/n)²]^(-1.5)
            where L(φ) = Lx cos²(φ) + Ly sin²(φ)
        """
        n_eff = max(int(n), 1)
        K = np.hypot(kx, ky)
        
        # Compute azimuth-dependent correlation length if anisotropic
        if self._is_anisotropic:
            # Compute azimuthal angle φ from spectral components
            if K > 1e-12:
                cos_phi = kx / K
                sin_phi = ky / K
                # L(φ) = Lx cos²(φ) + Ly sin²(φ)
                ell = self._lx * cos_phi**2 + self._ly * sin_phi**2
            else:
                # At K=0, use average correlation length
                ell = 0.5 * (self._lx + self._ly)
        else:
            # Isotropic case
            ell = float(correlation_length)
        
        # Compute spectrum: W^(n) = (L/n)² [1 + (KL/n)²]^(-1.5)
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
