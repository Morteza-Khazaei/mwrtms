"""High-level facade providing simplified APIs."""

from __future__ import annotations

from ..core.geometry import ScatteringGeometry
from ..core.polarization import PolarizationState
from ..core.wave import ElectromagneticWave
from ..factory import ScatteringModelFactory
from ..interface import ExponentialCorrelation, GaussianCorrelation, PowerLawCorrelation, SurfaceRoughness
from ..medium import IsotropicMedium
from ..scattering import ScatteringMechanism

__all__ = ["mwRTMsFacade"]


class mwRTMsFacade:
    """Facade exposing one-line helpers for common workflows."""

    @staticmethod
    def compute_soil_backscatter(
        model: str,
        frequency_ghz: float,
        incident_angle_deg: float,
        rms_height_cm: float,
        correlation_length_cm: float,
        soil_permittivity: complex,
        correlation: str = "exponential",
        polarization: PolarizationState = PolarizationState.VV,
        **model_kwargs,
    ) -> float:
        wave = ElectromagneticWave(frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=incident_angle_deg)
        corr_length_m = correlation_length_cm / 100.0
        rms_height_m = rms_height_cm / 100.0
        correlation_model = mwRTMsFacade._build_correlation(correlation, corr_length_m)
        roughness = SurfaceRoughness(rms_height_m, corr_length_m, correlation_model)
        air = IsotropicMedium(1.0, 290.0)
        soil = IsotropicMedium(soil_permittivity, 290.0)
        model_instance = ScatteringModelFactory.create(
            model,
            wave=wave,
            geometry=geometry,
            surface_roughness=roughness,
            **model_kwargs,
        )
        return model_instance.compute(air, soil, polarization)

    @staticmethod
    def build_model(model: str, **kwargs) -> ScatteringMechanism:
        return ScatteringModelFactory.create(model, **kwargs)

    @staticmethod
    def _build_correlation(kind: str, corr_length_m: float):
        kind_lower = kind.lower()
        if kind_lower == "gaussian":
            return GaussianCorrelation(corr_length_m)
        if kind_lower == "power":
            return PowerLawCorrelation(corr_length_m)
        return ExponentialCorrelation(corr_length_m)
