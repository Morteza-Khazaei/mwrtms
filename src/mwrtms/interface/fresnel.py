"""Fresnel reflection coefficients."""

from __future__ import annotations

import numpy as np

__all__ = ["FresnelCoefficients"]


class FresnelCoefficients:
    """Utility functions computing Fresnel reflection coefficients."""

    @staticmethod
    def reflection(eps1: complex, eps2: complex, theta: float, polarization: str) -> complex:
        cos_i = np.cos(theta)
        sin_i = np.sin(theta)

        sin_t = np.sqrt(eps1 / eps2) * sin_i
        cos_t = np.sqrt(1.0 - sin_t**2)

        if polarization.lower() == "h":
            return (cos_i - np.sqrt(eps2 / eps1) * cos_t) / (cos_i + np.sqrt(eps2 / eps1) * cos_t)
        if polarization.lower() == "v":
            return (eps2 * cos_i - np.sqrt(eps2 / eps1) * cos_t) / (eps2 * cos_i + np.sqrt(eps2 / eps1) * cos_t)
        raise ValueError("polarization must be 'h' or 'v'")
