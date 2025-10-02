"""High-level facade for mwRTMs."""

from __future__ import annotations

from typing import Dict

import numpy as np

from ..core import ElectromagneticWave, PolarizationState, ScatteringGeometry
from ..factory import ScatteringModelFactory
from ..interface import (
    AnisotropicRoughness,
    ExponentialCorrelation,
    GaussianCorrelation,
    IsotropicRoughness,
    PowerLawCorrelation,
)
from ..medium import Medium, SoilMedium

__all__ = ["mwRTMs"]


class mwRTMs:
    """Simplified facade exposing one-line backscatter computations."""

    @staticmethod
    def compute_soil_backscatter(
        model: str,
        frequency_ghz: float,
        incident_angle_deg: float,
        rms_height_cm: float,
        correlation_length_cm: float,
        correlation_type: str = "exponential",
        anisotropic: bool = False,
        correlation_length_y_cm: float | None = None,
        soil_moisture: float = 0.25,
        clay_fraction: float = 0.3,
        sand_fraction: float = 0.5,
        dielectric_model: str = "mironov",
        **model_kwargs,
    ) -> Dict[str, float]:
        """Return the backscatter in dB for HH, VV, HV channels."""

        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=incident_angle_deg, theta_s_deg=incident_angle_deg, phi_s_deg=180.0)

        if correlation_type == "exponential":
            corr_func = ExponentialCorrelation()
        elif correlation_type == "gaussian":
            corr_func = GaussianCorrelation()
        elif correlation_type == "powerlaw":
            corr_func = PowerLawCorrelation()
        else:
            raise ValueError(f"Unknown correlation type: {correlation_type}")

        if anisotropic:
            if correlation_length_y_cm is None:
                raise ValueError("correlation_length_y_cm required for anisotropic roughness")
            roughness = AnisotropicRoughness(
                rms_height_m=rms_height_cm / 100.0,
                correlation_length_x_m=correlation_length_cm / 100.0,
                correlation_length_y_m=correlation_length_y_cm / 100.0,
                correlation_function=corr_func,
            )
        else:
            roughness = IsotropicRoughness(
                rms_height_m=rms_height_cm / 100.0,
                correlation_length_m=correlation_length_cm / 100.0,
                correlation_function=corr_func,
            )

        soil = SoilMedium(
            moisture_m3m3=soil_moisture,
            clay_fraction=clay_fraction,
            sand_fraction=sand_fraction,
            dielectric_model=dielectric_model,
        )
        class _AirPlaceholder(Medium):
            def permittivity(self, frequency_hz: float) -> complex:  # pragma: no cover - trivial
                return 1.0 + 0.0j

        air = _AirPlaceholder()

        model_instance = ScatteringModelFactory.create(
            model,
            wave=wave,
            geometry=geometry,
            surface_roughness=roughness,
            **model_kwargs,
        )

        results: Dict[str, float] = {}
        for pol in (PolarizationState.HH, PolarizationState.VV, PolarizationState.HV):
            sigma = model_instance.compute(air, soil, pol)
            results[pol.value] = 10.0 * np.log10(sigma + 1e-30)
        return results
