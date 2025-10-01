"""Brightness temperature container."""

from __future__ import annotations

from ..core.polarization import StokesVector

__all__ = ["BrightnessTemperature"]


class BrightnessTemperature:
    """Encapsulated brightness container for passive models."""

    __slots__ = ("_tb_h", "_tb_v", "_frequency_hz", "_observation_angle_deg")
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, tb_h: float, tb_v: float, frequency_hz: float, observation_angle_deg: float) -> None:
        self._tb_h = float(tb_h)
        self._tb_v = float(tb_v)
        self._frequency_hz = float(frequency_hz)
        self._observation_angle_deg = float(observation_angle_deg)

    @property
    def TB_h(self) -> float:  # noqa: N802
        return self._tb_h

    @property
    def TB_v(self) -> float:  # noqa: N802
        return self._tb_v

    @property
    def frequency_hz(self) -> float:
        return self._frequency_hz

    @property
    def observation_angle_deg(self) -> float:
        return self._observation_angle_deg

    def polarization_ratio(self) -> float:
        return self._tb_v / self._tb_h

    def polarization_difference(self) -> float:
        return self._tb_v - self._tb_h

    def to_stokes(self) -> StokesVector:
        I = self._tb_h + self._tb_v
        Q = self._tb_v - self._tb_h
        return StokesVector(I=I, Q=Q, U=0.0, V=0.0)
