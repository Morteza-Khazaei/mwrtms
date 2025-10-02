"""Vegetation medium abstraction."""

from __future__ import annotations

from .dielectric import VegetationMaterialModel
from .base import Medium

__all__ = ["VegetationMedium"]


class VegetationMedium(Medium):
    """Homogeneous vegetation layer characterised by gravimetric moisture."""

    __slots__ = ("_gravimetric_moisture",)

    def __init__(self, gravimetric_moisture: float, temperature_k: float = 293.15) -> None:
        super().__init__(temperature_k)
        if not 0.0 <= gravimetric_moisture <= 5.0:
            raise ValueError("gravimetric_moisture must be within [0, 5]")
        self._gravimetric_moisture = float(gravimetric_moisture)
        self._dielectric_model = VegetationMaterialModel()

    @property
    def gravimetric_moisture(self) -> float:
        """Read-only access to the vegetation gravimetric moisture (g/g)."""

        return self._gravimetric_moisture

    def permittivity(self, frequency_hz: float) -> complex:
        """Return the complex permittivity (placeholder for Part 3)."""

        cached = self._dielectric_cache.get(frequency_hz)
        if cached is not None:
            return cached

        eps = self._dielectric_model.compute(
            frequency_hz=frequency_hz,
            gravimetric_moisture=self._gravimetric_moisture,
            temperature_k=self._temperature_k,
        )
        self._dielectric_cache[frequency_hz] = eps
        return eps
