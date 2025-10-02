"""Mironov soil dielectric model (MBSDM/GRMDM)."""

from __future__ import annotations

from typing import Final

import numpy as np

from .base import Medium

__all__ = ["MironovSoilMedium", "mironov_permittivity"]

_EPSILON_0: Final[float] = 8.854e-12  # F/m
_EINF: Final[float] = 4.9


def _spectral_parameters(clay_fraction: float) -> tuple[float, float, float, float, float, float, float, float, float, float]:
    c = clay_fraction
    nd = 1.634 - 0.539 * c + 0.2748 * c ** 2
    kd = 0.03952 - 0.04038 * c
    mvt = 0.02863 + 0.30673 * c
    e0b = 79.8 - 85.4 * c + 32.7 * c ** 2
    tb = 1.062e-11 + 3.450e-12 * c
    sb = 0.3112 + 0.467 * c
    su = 0.3631 + 1.217 * c
    e0u = 100.0
    tu = 8.5e-12
    return nd, kd, mvt, e0b, tb, sb, su, e0u, tu, c


def mironov_permittivity(moisture: float, clay_fraction: float, frequency_hz: float) -> complex:
    """Return the complex dielectric constant for mineral soils.

    Parameters
    ----------
    moisture:
        Volumetric soil moisture (0–1).
    clay_fraction:
        Clay content mass fraction (0–1).
    frequency_hz:
        Radar frequency in Hertz.

    Notes
    -----
    Implementation follows Mironov et al. (IEEE TGRS, 2009) using the
    mineralogically based spectroscopic dielectric model combined with the
    generalized refractive mixing approach.
    """

    if not 0.0 <= moisture <= 1.0:
        raise ValueError("moisture must be within [0, 1]")
    if not 0.0 <= clay_fraction <= 1.0:
        raise ValueError("clay_fraction must be within [0, 1]")
    if frequency_hz <= 0.0:
        raise ValueError("frequency_hz must be positive")

    omega = 2.0 * np.pi * frequency_hz
    nd, kd, mvt, e0b, tb, sb, su, e0u, tu, _ = _spectral_parameters(clay_fraction)

    eb = _EINF + (e0b - _EINF) / (1.0 - 1j * omega * tb) + 1j * sb / (omega * _EPSILON_0)
    eu = _EINF + (e0u - _EINF) / (1.0 - 1j * omega * tu) + 1j * su / (omega * _EPSILON_0)

    nb = np.sqrt(eb)
    nu = np.sqrt(eu)

    if moisture < mvt:
        nm = nd + (np.real(nb) - 1.0) * moisture
        km = kd + np.imag(nb) * moisture
    else:
        nm = nd + (np.real(nb) - 1.0) * mvt + (np.real(nu) - 1.0) * (moisture - mvt)
        km = kd + np.imag(nb) * mvt + np.imag(nu) * (moisture - mvt)

    emr = nm ** 2 - km ** 2
    emi = 2.0 * nm * km
    return complex(emr, emi)


class MironovSoilMedium(Medium):
    """Soil medium using the Mironov MBSDM dielectric formulation."""

    __slots__ = ("_moisture", "_clay_fraction", "_temperature")
    _mutable_slots: set[str] = set()

    def __init__(self, moisture: float, clay_fraction: float, *, temperature_k: float = 290.0) -> None:
        if not 0.0 <= moisture <= 1.0:
            raise ValueError("moisture must be within [0, 1]")
        if not 0.0 <= clay_fraction <= 1.0:
            raise ValueError("clay_fraction must be within [0, 1]")
        if temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive")
        self._moisture = float(moisture)
        self._clay_fraction = float(clay_fraction)
        self._temperature = float(temperature_k)

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    @property
    def moisture(self) -> float:
        return self._moisture

    @property
    def clay_fraction(self) -> float:
        return self._clay_fraction

    @property
    def temperature_k(self) -> float:
        return self._temperature

    def permittivity(self, frequency_hz: float, temperature_k: float | None = None) -> complex:
        return mironov_permittivity(self._moisture, self._clay_fraction, frequency_hz)

    def permeability(self, frequency_hz: float) -> complex:
        return complex(1.0)

    def physical_temperature(self) -> float:
        return self._temperature
