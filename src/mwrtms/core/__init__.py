"""Convenience imports for core mwRTMs abstractions."""

from .constants import EPSILON_0, MU_0, PLANCK_CONSTANT, SPEED_OF_LIGHT
from .geometry import ScatteringGeometry
from .polarization import (
    PolarizationState,
    default_polarization_order,
    normalize_polarization,
    normalize_polarization_sequence,
)
from .wave import ElectromagneticWave

from .radar_modes import (
    ObservationMode,
    RadarConfiguration,
    MonostaticConfiguration,
    BistaticConfiguration,
    RadarConfigurationFactory,
)

__all__ = [
    "SPEED_OF_LIGHT",
    "EPSILON_0",
    "MU_0",
    "PLANCK_CONSTANT",
    "ElectromagneticWave",
    "ScatteringGeometry",
    "PolarizationState",
    "default_polarization_order",
    "normalize_polarization",
    "normalize_polarization_sequence",
    "ObservationMode",
    "RadarConfiguration",
    "MonostaticConfiguration",
    "BistaticConfiguration",
    "RadarConfigurationFactory",
]
