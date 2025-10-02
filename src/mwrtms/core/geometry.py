"""Geometrical abstractions for active microwave scattering."""

from __future__ import annotations

from math import cos, radians, sin
from typing import Dict

import numpy as np

__all__ = ["ScatteringGeometry"]


class ScatteringGeometry:
    """Observation geometry encapsulating all angular descriptors.

    The class demonstrates encapsulation by keeping raw angles private and
    exposing read-only properties. Internally, frequently used conversions are
    cached to avoid repeated trigonometric evaluations.
    """

    __slots__ = ("_theta_i_deg", "_theta_s_deg", "_phi_i_deg", "_phi_s_deg", "_cache")

    def __init__(self, theta_i_deg: float, theta_s_deg: float | None = None, *, phi_s_deg: float = 180.0) -> None:
        self._validate_angle(theta_i_deg, "incident")
        if theta_s_deg is not None:
            self._validate_angle(theta_s_deg, "scattered")

        self._theta_i_deg = float(theta_i_deg)
        self._theta_s_deg = float(theta_s_deg) if theta_s_deg is not None else float(theta_i_deg)
        self._phi_i_deg = 0.0  # Incident azimuth fixed for backscatter setups
        self._phi_s_deg = float(phi_s_deg)
        self._cache: dict[str, float | Dict[str, np.ndarray]] = {}

    # ------------------------------------------------------------------
    # Encapsulated accessors
    # ------------------------------------------------------------------
    @property
    def theta_i_deg(self) -> float:
        """Incident zenith angle in degrees."""

        return self._theta_i_deg

    @property
    def theta_s_deg(self) -> float:
        """Scattered zenith angle in degrees."""

        return self._theta_s_deg

    @property
    def phi_s_deg(self) -> float:
        """Scattered azimuth angle in degrees."""

        return self._phi_s_deg

    @property
    def theta_i_rad(self) -> float:
        """Incident zenith angle in radians (cached)."""

        return self._cache.setdefault("theta_i_rad", radians(self._theta_i_deg))  # type: ignore[return-value]

    @property
    def theta_s_rad(self) -> float:
        """Scattered zenith angle in radians (cached)."""

        return self._cache.setdefault("theta_s_rad", radians(self._theta_s_deg))  # type: ignore[return-value]

    @property
    def phi_s_rad(self) -> float:
        """Scattered azimuth angle in radians (cached)."""

        return self._cache.setdefault("phi_s_rad", radians(self._phi_s_deg))  # type: ignore[return-value]

    @property
    def is_backscatter(self) -> bool:
        """Return ``True`` when the configuration is monostatic backscatter."""

        return np.isclose(self._theta_i_deg, self._theta_s_deg) and np.isclose((self._phi_s_deg % 360.0), 180.0)

    # ------------------------------------------------------------------
    # Derived quantities
    # ------------------------------------------------------------------
    def wave_vectors(self, wavenumber: float) -> Dict[str, np.ndarray]:
        """Return incident and scattered wave vectors for the given wavenumber.

        Parameters
        ----------
        wavenumber:
            Free-space wavenumber ``k`` in rad/m.
        """

        cache_key = f"wave_vectors_{float(wavenumber):.12g}"
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        ki = self._compute_wave_vector(wavenumber, self.theta_i_rad, 0.0)
        ks = self._compute_wave_vector(wavenumber, self.theta_s_rad, self.phi_s_rad)
        vectors = {"ki": ki, "ks": ks}
        self._cache[cache_key] = vectors
        return vectors

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _compute_wave_vector(k: float, theta: float, phi: float) -> np.ndarray:
        return k * np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])

    @staticmethod
    def _validate_angle(angle: float, name: str) -> None:
        if not 0.0 <= angle <= 90.0:
            raise ValueError(f"{name} angle must be within [0, 90] degrees")
