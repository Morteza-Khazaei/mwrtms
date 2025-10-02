"""Unit conversion helpers."""

from __future__ import annotations

import numpy as np

__all__ = ["cm_to_m", "ghz_to_hz", "linear_to_db"]


def cm_to_m(value_cm: float) -> float:
    return value_cm / 100.0


def ghz_to_hz(value_ghz: float) -> float:
    return value_ghz * 1e9


def linear_to_db(value_linear: float) -> float:
    return 10.0 * np.log10(value_linear + 1e-30)
