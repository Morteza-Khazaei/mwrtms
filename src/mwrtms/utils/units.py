"""Unit conversion helpers."""

from __future__ import annotations

__all__ = ["ghz_to_hz", "hz_to_ghz", "kelvin_to_celsius", "celsius_to_kelvin"]


def ghz_to_hz(value: float) -> float:
    return value * 1e9


def hz_to_ghz(value: float) -> float:
    return value / 1e9


def kelvin_to_celsius(value: float) -> float:
    return value - 273.15


def celsius_to_kelvin(value: float) -> float:
    return value + 273.15
