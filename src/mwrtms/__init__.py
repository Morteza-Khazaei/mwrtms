"""Foundation library for microwave radiative transfer models."""

from .__version__ import __version__
from .core.geometry import ScatteringGeometry
from .core.wave import ElectromagneticWave
from .core.polarization import PolarizationState
from .facade import mwRTMsFacade
from .factory import ScatteringModelFactory
from .interface import (
    CorrelationFunction,
    ExponentialCorrelation,
    GaussianCorrelation,
    PowerLawCorrelation,
    SurfaceRoughness,
)
from .medium import DielectricTensor, IsotropicMedium
from .scattering import (
    AIEMModel,
    Dubois95Model,
    I2EMModel,
    IEMModel,
    PRISM1Model,
    SMARTModel,
    SPM3DModel,
    SSRTModel,
)

__all__ = [
    "__version__",
    "ElectromagneticWave",
    "PolarizationState",
    "ScatteringGeometry",
    "SurfaceRoughness",
    "CorrelationFunction",
    "GaussianCorrelation",
    "ExponentialCorrelation",
    "PowerLawCorrelation",
    "DielectricTensor",
    "IsotropicMedium",
    "AIEMModel",
    "I2EMModel",
    "IEMModel",
    "PRISM1Model",
    "SPM3DModel",
    "SMARTModel",
    "Dubois95Model",
    "SSRTModel",
    "ScatteringModelFactory",
    "mwRTMsFacade",
]
