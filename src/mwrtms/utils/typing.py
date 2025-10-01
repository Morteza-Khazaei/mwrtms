"""Typing helpers."""

from __future__ import annotations

from typing import TypedDict

__all__ = ["SpectrumResult"]


class SpectrumResult(TypedDict):
    kx: float
    ky: float
    value: float
