"""Abstract scattering mechanism base classes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, Optional, TYPE_CHECKING, Union

from ..core.polarization import normalize_polarization
from ..result.scattering import ScatteringResult

__all__ = ["ScatteringMechanism"]

if TYPE_CHECKING:
    from ..core import PolarizationState
    from ..core.radar_modes import RadarConfiguration
    from ..medium import Medium


PolarizationInput = Optional[Union["PolarizationState", str, Iterable[Union["PolarizationState", str]]]]


class ScatteringMechanism(ABC):
    """Universal base class for scattering mechanisms."""

    MODEL_NAME = "ScatteringMechanism"

    def __init__(self, wave, geometry) -> None:
        self._wave = wave
        self._geometry = geometry

    @abstractmethod
    def compute(self, medium_above, medium_below, polarization) -> float:
        """Return the backscatter coefficient σ⁰ (linear power)."""

        raise NotImplementedError

    def run(
        self,
        medium_above: "Medium",
        medium_below: "Medium",
        polarizations: PolarizationInput = None,
        radar_config: Optional["RadarConfiguration"] = None,
    ) -> ScatteringResult:
        """Compute scattering for one or multiple polarizations.

        Parameters
        ----------
        medium_above, medium_below:
            Environmental media.
        polarizations:
            Optional polarization specifier (string, enum, iterable, or ``None`` for default order).
        radar_config:
            Optional configuration metadata to attach to the result.
        """

        pol_tuple = normalize_polarization(polarizations)
        data = {pol: self.compute(medium_above, medium_below, pol) for pol in pol_tuple}
        return ScatteringResult(
            data,
            getattr(self, "MODEL_NAME", self.__class__.__name__),
            self._wave,
            self._geometry,
            radar_config=radar_config,
        )

    def compute_with_config(
        self,
        medium_above: "Medium",
        medium_below: "Medium",
        polarization: PolarizationInput = None,
        config: Optional["RadarConfiguration"] = None,
    ) -> Union[ScatteringResult, float]:
        """Compute scattering using an explicit radar configuration."""

        if config is not None and not config.matches_geometry(self._geometry):
            raise ValueError("Radar configuration geometry does not match the scattering mechanism geometry")

        result = self.run(medium_above, medium_below, polarization, radar_config=config)
        pol_tuple = normalize_polarization(polarization)
        if len(pol_tuple) == 1:
            return result[pol_tuple[0]]
        return result

    @property
    def wave(self):
        return self._wave

    @property
    def geometry(self):
        return self._geometry
