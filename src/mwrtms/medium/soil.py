"""Soil medium abstraction."""

from __future__ import annotations

from ..dielectric import DobsonModel, MironovModel
from .base import Medium

__all__ = ["SoilMedium"]


class SoilMedium(Medium):
    """Homogeneous soil medium with encapsulated texture parameters."""

    __slots__ = (
        "_moisture",
        "_clay",
        "_sand",
        "_silt",
        "_dielectric_model",
    )

    def __init__(
        self,
        moisture_m3m3: float,
        clay_fraction: float,
        sand_fraction: float,
        temperature_k: float = 293.15,
        dielectric_model: str = "mironov",
    ) -> None:
        super().__init__(temperature_k)
        self._validate_parameters(moisture_m3m3, clay_fraction, sand_fraction)
        self._moisture = float(moisture_m3m3)
        self._clay = float(clay_fraction)
        self._sand = float(sand_fraction)
        self._silt = 1.0 - self._clay - self._sand

        model_name = dielectric_model.lower()
        if model_name == "dobson":
            self._dielectric_model = DobsonModel()
        elif model_name == "mironov":
            self._dielectric_model = MironovModel()
        else:
            raise ValueError(f"Unknown dielectric model: {dielectric_model}")

    # ------------------------------------------------------------------
    # Encapsulated read-only properties
    # ------------------------------------------------------------------
    @property
    def moisture(self) -> float:
        """Volumetric soil moisture (m³/m³)."""

        return self._moisture

    @property
    def clay_fraction(self) -> float:  # type: ignore[override]
        """Clay fraction of the soil texture (0–1)."""

        return self._clay

    @property
    def sand_fraction(self) -> float:  # type: ignore[override]
        """Sand fraction of the soil texture (0–1)."""

        return self._sand

    @property
    def silt_fraction(self) -> float:
        """Silt fraction of the soil texture (0–1)."""

        return self._silt

    @property
    def dielectric_model_name(self) -> str:
        """Name of the dielectric model selected for this soil."""

        return self._dielectric_model.__class__.__name__.lower()

    # ------------------------------------------------------------------
    # Medium interface
    # ------------------------------------------------------------------
    def permittivity(self, frequency_hz: float) -> complex:
        """Return the complex permittivity (placeholder for Part 3)."""

        cached = self._dielectric_cache.get(frequency_hz)
        if cached is not None:
            return cached

        eps = self._dielectric_model.compute(
            frequency_hz=frequency_hz,
            moisture=self._moisture,
            clay=self._clay,
            sand=self._sand,
            temperature_k=self._temperature_k,
        )
        self._dielectric_cache[frequency_hz] = eps
        return eps

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _validate_parameters(moisture: float, clay: float, sand: float) -> None:
        if not 0.0 <= moisture <= 1.0:
            raise ValueError("moisture must be within [0, 1]")
        if not 0.0 <= clay <= 1.0:
            raise ValueError("clay_fraction must be within [0, 1]")
        if not 0.0 <= sand <= 1.0:
            raise ValueError("sand_fraction must be within [0, 1]")
        if clay + sand > 1.0:
            raise ValueError("clay_fraction + sand_fraction must not exceed 1.0")
