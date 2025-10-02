"""Reusable scattering scenario helper."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Union

from ..core import ElectromagneticWave, PolarizationState, RadarConfiguration
from ..core.polarization import normalize_polarization
from ..medium import HomogeneousMedium, Medium
from ..medium.surface import Surface, build_surface_from_statistics
from ..factory import ScatteringModelFactory
from ..result.scattering import ScatteringResult


PolarizationSpec = Optional[Union[PolarizationState, str, Iterable[Union[PolarizationState, str]]]]


@dataclass
class ScatteringScene:
    """Encapsulate the ingredients required for a scattering simulation."""

    radar_config: RadarConfiguration
    wave: ElectromagneticWave
    medium_above: Medium
    medium_below: Medium
    surface: Optional[Surface] = None

    @classmethod
    def soil_scene(
        cls,
        radar_config: RadarConfiguration,
        frequency_ghz: float,
        soil_permittivity: complex,
        *,
        rms_height_cm: Optional[float] = None,
        correlation_length_cm: Optional[float] = None,
        correlation_length_y_cm: Optional[float] = None,
        correlation: str = "exponential",
        air_permittivity: complex = 1.0 + 0.0j,
    ) -> "ScatteringScene":
        """Create a soil scene with homogeneous media and optional surface statistics."""

        corr_type = correlation.lower()
        if corr_type not in {"exponential", "gaussian", "powerlaw"}:
            raise ValueError(f"Unknown correlation type: {correlation}")

        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        air = HomogeneousMedium(air_permittivity)
        soil = HomogeneousMedium(soil_permittivity)

        surface: Optional[Surface] = None
        if rms_height_cm is not None and correlation_length_cm is not None:
            surface = build_surface_from_statistics(
                rms_height_cm / 100.0,
                correlation_length_cm / 100.0,
                correlation_length_y_m=correlation_length_y_cm / 100.0 if correlation_length_y_cm else None,
                correlation_type=corr_type,
            )

        return cls(radar_config, wave, air, soil, surface)

    def with_surface(self, surface: Surface) -> "ScatteringScene":
        """Return a copy of the scene with an explicit surface assigned."""

        return ScatteringScene(
            radar_config=self.radar_config,
            wave=self.wave,
            medium_above=self.medium_above,
            medium_below=self.medium_below,
            surface=surface,
        )

    def ensure_surface(
        self,
        surface: Optional[Surface] = None,
        *,
        rms_height_m: Optional[float] = None,
        correlation_length_m: Optional[float] = None,
        correlation_length_y_m: Optional[float] = None,
        correlation_type: str = "exponential",
    ) -> Surface:
        """Return a surface either from arguments or stored on the scene."""

        if surface is not None:
            return surface
        if self.surface is not None:
            return self.surface
        if rms_height_m is not None and correlation_length_m is not None:
            return build_surface_from_statistics(
                rms_height_m,
                correlation_length_m,
                correlation_length_y_m=correlation_length_y_m,
                correlation_type=correlation_type,
            )
        raise ValueError("A surface must be supplied or constructed via statistics")

    def run_model(
        self,
        model: str,
        *,
        polarizations: PolarizationSpec = None,
        surface: Optional[Surface] = None,
        model_kwargs: Optional[dict] = None,
        **extra_kwargs,
    ) -> ScatteringResult:
        """Create the model and execute it for the requested polarizations."""

        surface_obj = self.ensure_surface(surface, **extra_kwargs)
        kwargs = dict(model_kwargs or {})
        model_instance = ScatteringModelFactory.create_with_radar_config(
            model,
            config=self.radar_config,
            wave=self.wave,
            surface=surface_obj,
            **kwargs,
        )
        return model_instance.run(
            self.medium_above,
            self.medium_below,
            polarizations=polarizations,
            radar_config=self.radar_config,
        )

    def run_scalar(
        self,
        model: str,
        polarization: Union[PolarizationState, str],
        **model_kwargs,
    ) -> float:
        """Convenience helper returning a single polarization backscatter."""

        result = self.run_model(model, polarizations=polarization, model_kwargs=model_kwargs)
        pol = normalize_polarization(polarization)[0]
        return result[pol]
