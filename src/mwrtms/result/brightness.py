"""Brightness temperature result container."""

from __future__ import annotations

from ..core.geometry import ScatteringGeometry
from ..core.wave import ElectromagneticWave

__all__ = ["BrightnessTemperatureResult"]


class BrightnessTemperatureResult:
    """Encapsulated brightness temperature container."""

    __slots__ = (
        "_tb_h",
        "_tb_v",
        "_emissivity_h",
        "_emissivity_v",
        "_physical_temperature",
        "_model_name",
        "_wave",
        "_geometry",
    )
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(
        self,
        tb_h: float,
        tb_v: float,
        emissivity_h: float,
        emissivity_v: float,
        physical_temperature: float,
        model_name: str,
        wave: ElectromagneticWave,
        geometry: ScatteringGeometry,
    ) -> None:
        self._tb_h = float(tb_h)
        self._tb_v = float(tb_v)
        self._emissivity_h = float(emissivity_h)
        self._emissivity_v = float(emissivity_v)
        self._physical_temperature = float(physical_temperature)
        self._model_name = model_name
        self._wave = wave
        self._geometry = geometry

    @property
    def TB_h(self) -> float:  # noqa: N802
        return self._tb_h

    @property
    def TB_v(self) -> float:  # noqa: N802
        return self._tb_v

    @property
    def emissivity_h(self) -> float:
        return self._emissivity_h

    @property
    def emissivity_v(self) -> float:
        return self._emissivity_v

    @property
    def physical_temperature(self) -> float:
        return self._physical_temperature

    @property
    def model_name(self) -> str:
        return self._model_name

    @property
    def wave(self) -> ElectromagneticWave:
        return self._wave

    @property
    def geometry(self) -> ScatteringGeometry:
        return self._geometry

    def polarization_ratio(self) -> float:
        return self._tb_v / self._tb_h

    def polarization_difference(self) -> float:
        return self._tb_v - self._tb_h
