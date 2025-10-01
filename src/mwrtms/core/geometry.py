"""Geometrical abstractions for scattering and emission."""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np

__all__ = ["ScatteringGeometry"]


class ScatteringGeometry:
    """Observation geometry for scattering/emission problems."""

    __slots__ = (
        "_theta_i_deg",
        "_theta_s_deg",
        "_phi_i_deg",
        "_phi_s_deg",
        "_theta_i_rad_cache",
        "_theta_s_rad_cache",
        "_wave_vector_cache",
    )
    _mutable_slots = {"_theta_i_rad_cache", "_theta_s_rad_cache", "_wave_vector_cache"}

    def __setattr__(self, name, value):
        if name in self.__slots__ and name not in self._mutable_slots and hasattr(self, name):
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(
        self,
        theta_i_deg: float,
        theta_s_deg: Optional[float] = None,
        phi_i_deg: float = 0.0,
        phi_s_deg: float = 180.0,
    ) -> None:
        self._validate_angle(theta_i_deg, "incident")
        if theta_s_deg is not None:
            self._validate_angle(theta_s_deg, "scattered")
        self._theta_i_deg = float(theta_i_deg)
        self._theta_s_deg = float(theta_s_deg) if theta_s_deg is not None else float(theta_i_deg)
        self._phi_i_deg = float(phi_i_deg)
        self._phi_s_deg = float(phi_s_deg)
        self._theta_i_rad_cache: Optional[float] = None
        self._theta_s_rad_cache: Optional[float] = None
        self._wave_vector_cache: Optional[tuple[float, Dict[str, np.ndarray]]] = None

    @property
    def theta_i_deg(self) -> float:
        return self._theta_i_deg

    @property
    def theta_s_deg(self) -> float:
        return self._theta_s_deg

    @property
    def phi_i_deg(self) -> float:
        return self._phi_i_deg

    @property
    def phi_s_deg(self) -> float:
        return self._phi_s_deg

    @property
    def theta_i(self) -> float:
        if self._theta_i_rad_cache is None:
            self._theta_i_rad_cache = np.deg2rad(self._theta_i_deg)
        return self._theta_i_rad_cache

    @property
    def theta_s(self) -> float:
        if self._theta_s_rad_cache is None:
            self._theta_s_rad_cache = np.deg2rad(self._theta_s_deg)
        return self._theta_s_rad_cache

    @property
    def phi_i(self) -> float:
        return np.deg2rad(self._phi_i_deg)

    @property
    def phi_s(self) -> float:
        return np.deg2rad(self._phi_s_deg)

    @property
    def is_backscatter(self) -> bool:
        return bool(
            np.isclose(self._theta_i_deg, self._theta_s_deg)
            and np.isclose(abs(self._phi_s_deg - self._phi_i_deg), 180.0)
        )

    @property
    def is_nadir(self) -> bool:
        return bool(np.isclose(self._theta_i_deg, 0.0))

    def wave_vectors(self, wavenumber: float) -> Dict[str, np.ndarray]:
        cache_key = float(wavenumber)
        if self._wave_vector_cache is None or self._wave_vector_cache[0] != cache_key:
            ki = self._compute_wave_vector(wavenumber, self.theta_i, self._phi_i_deg)
            ks = self._compute_wave_vector(wavenumber, self.theta_s, self._phi_s_deg)
            self._wave_vector_cache = (cache_key, {"ki": ki, "ks": ks})
        return self._wave_vector_cache[1]

    def _compute_wave_vector(self, k: float, theta: float, phi_deg: float) -> np.ndarray:
        phi = np.deg2rad(phi_deg)
        return k * np.array([
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ])

    def _validate_angle(self, angle: float, name: str) -> None:
        if not 0.0 <= angle <= 90.0:
            raise ValueError(f"{name} angle must be within [0, 90] degrees")
