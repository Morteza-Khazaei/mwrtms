"""Fresnel reflection and transmission coefficients."""

from __future__ import annotations

import cmath

import numpy as np

__all__ = ["FresnelCoefficients"]


class FresnelCoefficients:
    """Utility functions for Fresnel reflection, transmission, emissivity."""

    @staticmethod
    def _transmitted_angle(eps1: complex, eps2: complex, theta_i: float) -> float:
        sin_t_sq = (eps1 / eps2) * np.sin(theta_i) ** 2
        return float(np.arcsin(np.sqrt(sin_t_sq)))

    @staticmethod
    def reflection(eps1: complex, eps2: complex, theta_i: float, polarization: str) -> complex:
        pol = polarization.lower()
        cos_i = np.cos(theta_i)
        sin_i = np.sin(theta_i)
        eta = cmath.sqrt(eps2 / eps1)
        sin_t = sin_i / eta
        cos_t = cmath.sqrt(1.0 - sin_t**2)
        if pol == "h":
            return (cos_i - eta * cos_t) / (cos_i + eta * cos_t)
        if pol == "v":
            return (eps2 * cos_i - eps1 * cos_t) / (eps2 * cos_i + eps1 * cos_t)
        raise ValueError("polarization must be 'h' or 'v'")

    @staticmethod
    def transmission(eps1: complex, eps2: complex, theta_i: float, polarization: str) -> complex:
        return 1.0 + FresnelCoefficients.reflection(eps1, eps2, theta_i, polarization)

    @staticmethod
    def emissivity(eps1: complex, eps2: complex, theta_i: float, polarization: str) -> float:
        refl = FresnelCoefficients.reflection(eps1, eps2, theta_i, polarization)
        return 1.0 - abs(refl) ** 2
