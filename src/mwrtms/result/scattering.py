"""Scattering result container."""

from __future__ import annotations

from typing import Dict, Mapping, Optional

import numpy as np

from ..core.geometry import ScatteringGeometry
from ..core.wave import ElectromagneticWave

__all__ = ["ScatteringResult"]


class ScatteringResult:
    """Immutable scattering result that encapsulates model outputs."""

    __slots__ = ("_data", "_model_name", "_wave", "_geometry", "_components")
    _mutable_slots = set()

    def __init__(
        self,
        data: Mapping[str, float],
        model_name: str,
        wave: ElectromagneticWave,
        geometry: ScatteringGeometry,
        components: Optional[Mapping[str, float]] = None,
    ) -> None:
        normalized = {pol.lower(): float(value) for pol, value in data.items()}
        for key, value in normalized.items():
            if value < 0:
                raise ValueError(f"Scattering coefficient for {key} must be non-negative")
        self._data = normalized
        self._model_name = model_name
        self._wave = wave
        self._geometry = geometry
        self._components = dict(components) if components is not None else None

    @property
    def model_name(self) -> str:
        return self._model_name

    @property
    def wave(self) -> ElectromagneticWave:
        return self._wave

    @property
    def geometry(self) -> ScatteringGeometry:
        return self._geometry

    @property
    def components(self) -> Optional[Mapping[str, float]]:
        return None if self._components is None else dict(self._components)

    def __getitem__(self, pol: str) -> float:
        return self._data[pol.lower()]

    def to_db(self) -> "ScatteringResult":
        converted: Dict[str, float] = {pol: 10.0 * np.log10(val + 1e-30) for pol, val in self._data.items()}
        return ScatteringResult(converted, self._model_name, self._wave, self._geometry, self._components)

    def polarization_ratio(self, numerator: str, denominator: str) -> float:
        return 10.0 * np.log10(self[numerator] / self[denominator])
