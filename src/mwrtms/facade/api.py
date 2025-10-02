"""High-level facade for mwRTMs."""

from __future__ import annotations

from typing import Iterable, Union

from ..core import PolarizationState, RadarConfiguration
from ..core.polarization import normalize_polarization
from ..integration import ScatteringScene
from ..result.scattering import ScatteringResult

__all__ = ["mwRTMs"]


PolarizationInput = Union[PolarizationState, str, Iterable[Union[PolarizationState, str]], None]


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
        polarization: PolarizationInput = None,
        **model_kwargs,
    ) -> Union[float, ScatteringResult]:
        """Run a surface scattering model backed by the :class:`ScatteringScene` helper.

        When ``polarization`` is ``None`` all default channels are simulated and a
        :class:`ScatteringResult` is returned. Passing a single polarization keeps
        the historical scalar return behaviour.
        """

        scene = ScatteringScene.soil_scene(
            radar_config=radar_config,
            frequency_ghz=frequency_ghz,
            soil_permittivity=soil_permittivity,
            rms_height_cm=rms_height_cm,
            correlation_length_cm=correlation_length_cm,
            correlation=correlation,
        )

        result = scene.run_model(model, polarizations=polarization, model_kwargs=model_kwargs)
        pol_tuple = normalize_polarization(polarization)
        if len(pol_tuple) == 1:
            return result[pol_tuple[0]]
        return result
