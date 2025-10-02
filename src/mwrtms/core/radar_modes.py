"""Radar observation mode abstractions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, List, Sequence

import numpy as np

from .geometry import ScatteringGeometry

__all__ = [
    "ObservationMode",
    "RadarConfiguration",
    "MonostaticConfiguration",
    "BistaticConfiguration",
    "RadarConfigurationFactory",
]


class ObservationMode(str, Enum):
    """Supported radar observation modes."""

    MONOSTATIC = "monostatic"
    BISTATIC = "bistatic"


@dataclass(frozen=True)
class RadarConfiguration:
    """Base radar configuration identified by a geometry and observation mode."""

    geometry: ScatteringGeometry
    mode: ObservationMode
    description: str | None = None

    def matches_geometry(self, geometry: ScatteringGeometry, *, atol: float = 1e-6) -> bool:
        """Return ``True`` when the supplied geometry matches this configuration."""

        return (
            np.isclose(self.geometry.theta_i_deg, geometry.theta_i_deg, atol=atol)
            and np.isclose(self.geometry.theta_s_deg, geometry.theta_s_deg, atol=atol)
            and np.isclose((self.geometry.phi_s_deg % 360.0), (geometry.phi_s_deg % 360.0), atol=atol)
        )

    def bistatic_angle(self) -> float:
        """Return the bistatic angle in degrees if defined."""

        if self.mode != ObservationMode.BISTATIC:
            raise ValueError("Bistatic angle only defined for bistatic configurations")
        return abs(self.geometry.theta_s_deg - self.geometry.theta_i_deg)

    @property
    def is_monostatic(self) -> bool:
        return self.mode == ObservationMode.MONOSTATIC

    @property
    def is_bistatic(self) -> bool:
        return self.mode == ObservationMode.BISTATIC


class MonostaticConfiguration(RadarConfiguration):
    """Monostatic backscatter configuration."""

    def __init__(self, theta_deg: float, *, phi_s_deg: float = 180.0, description: str | None = None) -> None:
        geometry = ScatteringGeometry(theta_i_deg=theta_deg, theta_s_deg=theta_deg, phi_s_deg=phi_s_deg)
        super().__init__(geometry=geometry, mode=ObservationMode.MONOSTATIC, description=description)

    @property
    def look_angle_deg(self) -> float:
        return self.geometry.theta_i_deg


class BistaticConfiguration(RadarConfiguration):
    """Bistatic configuration with arbitrary incident/scattered directions."""

    def __init__(
        self,
        theta_i_deg: float,
        theta_s_deg: float,
        *,
        phi_s_deg: float = 0.0,
        description: str | None = None,
    ) -> None:
        geometry = ScatteringGeometry(theta_i_deg=theta_i_deg, theta_s_deg=theta_s_deg, phi_s_deg=phi_s_deg)
        super().__init__(geometry=geometry, mode=ObservationMode.BISTATIC, description=description)


class RadarConfigurationFactory:
    """Factory helpers for radar configuration creation."""

    @staticmethod
    def create_monostatic(
        theta_deg: float,
        *,
        phi_s_deg: float = 180.0,
        description: str | None = None,
    ) -> MonostaticConfiguration:
        return MonostaticConfiguration(theta_deg, phi_s_deg=phi_s_deg, description=description)

    @staticmethod
    def create_bistatic(
        theta_i_deg: float,
        theta_s_deg: float,
        *,
        phi_s_deg: float = 0.0,
        description: str | None = None,
    ) -> BistaticConfiguration:
        return BistaticConfiguration(
            theta_i_deg=theta_i_deg,
            theta_s_deg=theta_s_deg,
            phi_s_deg=phi_s_deg,
            description=description,
        )

    @staticmethod
    def create_multi_angle_monostatic(
        angles_deg: Sequence[float] | Iterable[float],
        *,
        phi_s_deg: float = 180.0,
        description: str | None = None,
    ) -> List[MonostaticConfiguration]:
        base_desc = description or "Monostatic"
        return [
            MonostaticConfiguration(angle, phi_s_deg=phi_s_deg, description=f"{base_desc} θ={angle:.1f}°")
            for angle in angles_deg
        ]

    @staticmethod
    def from_geometry(
        geometry: ScatteringGeometry,
        *,
        mode: ObservationMode | None = None,
        description: str | None = None,
    ) -> RadarConfiguration:
        detected_mode = mode or (ObservationMode.MONOSTATIC if geometry.is_backscatter else ObservationMode.BISTATIC)
        if detected_mode == ObservationMode.MONOSTATIC:
            return MonostaticConfiguration(geometry.theta_i_deg, phi_s_deg=geometry.phi_s_deg, description=description)
        return BistaticConfiguration(
            geometry.theta_i_deg,
            geometry.theta_s_deg,
            phi_s_deg=geometry.phi_s_deg,
            description=description,
        )
