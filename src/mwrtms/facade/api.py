"""High-level facade for mwRTMs."""

from __future__ import annotations

from ..core import ElectromagneticWave, PolarizationState, RadarConfiguration
from ..factory import ScatteringModelFactory
from ..medium.surface import build_surface_from_statistics
from ..medium import Medium

__all__ = ["mwRTMs"]


class mwRTMs:
    """Simplified facade exposing one-line backscatter computations."""

    @staticmethod
    def compute_soil_backscatter(
        model: str,
        radar_config: RadarConfiguration,
        frequency_ghz: float,
        rms_height_cm: float,
        correlation_length_cm: float,
        soil_permittivity: complex,
        correlation: str = "exponential",
        polarization: PolarizationState = PolarizationState.VV,
        **model_kwargs,
    ) -> float:
        """Return the backscatter coefficient (linear) for a radar configuration."""

        corr_type = correlation.lower()
        if corr_type not in {"exponential", "gaussian", "powerlaw"}:
            raise ValueError(f"Unknown correlation type: {correlation}")

        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        surface = build_surface_from_statistics(
            rms_height_cm / 100.0,
            correlation_length_cm / 100.0,
            correlation_type=corr_type,
        )

        class _ConstantMedium(Medium):
            def __init__(self, permittivity: complex) -> None:
                super().__init__(temperature_k=293.15)
                self._eps = permittivity

            def permittivity(self, frequency_hz: float) -> complex:  # pragma: no cover - trivial
                return self._eps

        air = _ConstantMedium(1.0 + 0.0j)
        soil = _ConstantMedium(soil_permittivity)

        model_instance = ScatteringModelFactory.create_with_radar_config(
            model,
            config=radar_config,
            wave=wave,
            surface=surface,
            **model_kwargs,
        )

        return model_instance.compute_with_config(air, soil, polarization, radar_config)
