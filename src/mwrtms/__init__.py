"""Microwave Radiative Transfer Models foundation package."""

from .__version__ import __version__
from .core.geometry import ScatteringGeometry
from .core.wave import ElectromagneticWave
from .core.polarization import PolarizationState
from .medium import (
    Medium,
    SoilMedium,
    VegetationMedium,
    MironovSoilMedium,
    mironov_permittivity,
)
from .interface import (
    ExponentialCorrelation,
    GaussianCorrelation,
    PowerLawCorrelation,
    IsotropicRoughness,
    AnisotropicRoughness,
)
from .scattering import (
    ScatteringMechanism,
    SurfaceScattering,
    AIEMModel,
    SPMModel,
    VolumeScattering,
    SSRTModel,
    CanopyProperties,
)
from .factory import ScatteringModelFactory
from .facade import mwRTMs

__all__ = [
    "__version__",
    "ElectromagneticWave",
    "ScatteringGeometry",
    "PolarizationState",
    "Medium",
    "SoilMedium",
    "VegetationMedium",
    "MironovSoilMedium",
    "mironov_permittivity",
    "ExponentialCorrelation",
    "GaussianCorrelation",
    "PowerLawCorrelation",
    "IsotropicRoughness",
    "AnisotropicRoughness",
    "ScatteringMechanism",
    "SurfaceScattering",
    "AIEMModel",
    "SPMModel",
    "VolumeScattering",
    "SSRTModel",
    "CanopyProperties",
    "ScatteringModelFactory",
    "mwRTMs",
]
