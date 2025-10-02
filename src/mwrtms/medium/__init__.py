"""Medium hierarchy exports."""

from .base import Medium
from .soil import SoilMedium
from .vegetation import VegetationMedium
from .dielectric import DobsonModel, MironovModel, VegetationMaterialModel
from .homogeneous import HomogeneousMedium
from .interface import FresnelCoefficients
from .surface import (
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

__all__ = [
    "Medium",
    "SoilMedium",
    "VegetationMedium",
    "DobsonModel",
    "MironovModel",
    "VegetationMaterialModel",
    "HomogeneousMedium",
    "Surface",
    "SurfaceAnalyzer",
    "SurfaceGenerator",
    "SyntheticSurfaceGenerator",
    "MeasuredSurfaceLoader",
    "CompositeSurfaceGenerator",
    "build_surface_from_statistics",
    "CorrelationFunction",
    "Exponential",
    "Gaussian",
    "PowerLaw",
    "FresnelCoefficients",
    "MironovSoilMedium",
    "mironov_permittivity",
]


class MironovSoilMedium(SoilMedium):
    """Convenience wrapper that locks SoilMedium to the Mironov dielectric model."""

    def __init__(
        self,
        moisture_m3m3: float,
        clay_fraction: float,
        sand_fraction: float,
        temperature_k: float = 293.15,
    ) -> None:
        super().__init__(
            moisture_m3m3=moisture_m3m3,
            clay_fraction=clay_fraction,
            sand_fraction=sand_fraction,
            temperature_k=temperature_k,
            dielectric_model="mironov",
        )


def mironov_permittivity(
    frequency_hz: float,
    *,
    moisture: float,
    clay: float,
    sand: float,
    temperature_k: float = 293.15,
) -> complex:
    """Return Mironov et al. (2009) soil permittivity for convenience."""

    model = MironovModel()
    return model.compute(
        frequency_hz=frequency_hz,
        moisture=moisture,
        clay=clay,
        sand=sand,
        temperature_k=temperature_k,
    )
