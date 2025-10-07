"""Correlation function helpers for IEM-family models."""

from __future__ import annotations

import math
from typing import Iterable, Tuple

import numpy as np
from scipy.special import gamma, kv

__all__ = [
    "single_scale_spectrum_weights",
    "multiscale_gaussian_weights",
    "rms_slope_from_correlation",
    "auto_spectral_order",
]

_CORRELATION_ALIASES = {
    1: "exponential",
    2: "gaussian",
    3: "powerlaw",
    4: "x-exponential",
    "exp": "exponential",
    "gauss": "gaussian",
    "power": "powerlaw",
    "powerlaw": "powerlaw",
    "xpower": "powerlaw",
    "x-exponential": "x-exponential",
}


def _normalize_correlation(correlation: str | int) -> str:
    label = correlation
    if isinstance(label, str):
        label = label.lower()
    return _CORRELATION_ALIASES.get(label, str(label))


def rms_slope_from_correlation(
    correlation: str | int,
    sigma: float,
    corr_length: float,
    power_exponent: float = 1.5,
) -> float:
    """Return RMS slope (dimensionless) for the requested correlation."""

    corr = _normalize_correlation(correlation)
    ratio = sigma / corr_length
    if corr == "exponential" or corr == "x-exponential":
        return ratio
    if corr == "gaussian":
        return math.sqrt(2.0) * ratio
    if corr == "powerlaw":
        return math.sqrt(2.0 * power_exponent) * ratio
    raise NotImplementedError(f"Unsupported correlation type '{correlation}'")


def single_scale_spectrum_weights(
    wvnb: float,
    sigma: float,
    corr_length: float,
    n_terms: int,
    correlation: str | int,
    power_exponent: float = 1.5,
) -> Tuple[np.ndarray, float]:
    """Compute the I2EM roughness spectrum weights for a single-scale surface."""

    corr = _normalize_correlation(correlation)
    w = np.zeros(n_terms, dtype=float)
    L = corr_length
    x = wvnb * L

    if corr == "exponential":
        n = np.arange(1, n_terms + 1)
        w = (L ** 2) / (n ** 2) * (1.0 + (x / n) ** 2) ** (-1.5)
    elif corr == "gaussian":
        n = np.arange(1, n_terms + 1)
        w = (L ** 2) / (2.0 * n) * np.exp(-(x ** 2) / (4.0 * n))
    elif corr == "powerlaw":
        n = np.arange(1, n_terms + 1)
        if np.isclose(x, 0.0):
            w = (L ** 2) / (3.0 * n - 2.0)
        else:
            nu = power_exponent * n
            w = (
                (L ** 2)
                * (x ** (-1.0 + nu))
                * kv(1.0 - nu, x)
                / (2.0 ** (nu - 1.0) * gamma(nu))
            )
    else:
        raise NotImplementedError(f"Correlation '{correlation}' not yet supported")

    slope = rms_slope_from_correlation(corr, sigma, corr_length, power_exponent)
    return w, slope


def multiscale_gaussian_weights(
    sigmas: Iterable[float],
    lengths: Iterable[float],
    wavenumber: float,
    n_terms: int,
) -> np.ndarray:
    """Return Gaussian spectrum weights for multi-scale I2EM surfaces."""

    sigmas = np.asarray(list(sigmas), dtype=float)
    lengths = np.asarray(list(lengths), dtype=float)
    if sigmas.shape != lengths.shape:
        raise ValueError("sigmas and lengths must have matching shapes")

    sigma_sq = sigmas ** 2
    total_sigma_sq = np.sum(sigma_sq)
    if total_sigma_sq <= 0:
        raise ValueError("Total variance must be positive")

    n = np.arange(1, n_terms + 1)
    weights = np.zeros(n_terms, dtype=float)
    for sig_sq, L in zip(sigma_sq, lengths):
        L_m = L
        x = wavenumber * L_m
        weights += (sig_sq / total_sigma_sq) * (L_m ** 2) / (2.0 * n) * np.exp(-(x ** 2) / (4.0 * n))
    return weights


def auto_spectral_order(ks: float, cos_theta_i: float, cos_theta_s: float, tolerance: float = 1e-8) -> int:
    """Select the number of spectral terms required for convergence."""

    order = 1
    error = 1.0e8
    argument = (ks ** 2) * (cos_theta_i + cos_theta_s) ** 2
    while error > tolerance:
        order += 1
        error = (argument ** order) / math.factorial(order)
    return order
