"""Microwave Radiative Transfer Models foundation package."""

from .__version__ import __version__
from .core.geometry import ScatteringGeometry
from .core.wave import ElectromagneticWave
from .core.polarization import PolarizationState

from .core.radar_modes import (
    ObservationMode,
    RadarConfiguration,
    MonostaticConfiguration,
    BistaticConfiguration,
    RadarConfigurationFactory,
)
from .medium import (
    Medium,
    SoilMedium,
    VegetationMedium,
    MironovSoilMedium,
    mironov_permittivity,
    FresnelCoefficients,
    HomogeneousMedium,
)
from .medium.surface import (
    Surface,
    SurfaceAnalyzer,
    SurfaceGenerator,
    SyntheticSurfaceGenerator,
    MeasuredSurfaceLoader,
    CompositeSurfaceGenerator,
    build_surface_from_statistics,
    CorrelationFunction,
    Exponential,
    Gaussian,
    PowerLaw,
)
from .integration import ScatteringScene
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
from .facade import mwRTMs, mwRTMsFacade

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
    "HomogeneousMedium",
    "ScatteringScene",
    "FresnelCoefficients",
    "ObservationMode",
    "RadarConfiguration",
    "MonostaticConfiguration",
    "BistaticConfiguration",
    "RadarConfigurationFactory",
    "CorrelationFunction",
    "Exponential",
    "Gaussian",
    "PowerLaw",
    "Surface",
    "SurfaceAnalyzer",
    "SurfaceGenerator",
    "SyntheticSurfaceGenerator",
    "MeasuredSurfaceLoader",
    "CompositeSurfaceGenerator",
    "build_surface_from_statistics",
    "ScatteringScene",
    "ScatteringMechanism",
    "SurfaceScattering",
    "AIEMModel",
    "SPMModel",
    "VolumeScattering",
    "SSRTModel",
    "CanopyProperties",
    "ScatteringModelFactory",
    "mwRTMs",
    "mwRTMsFacade",
]
