"""Placeholder SMART model while implementation is unavailable."""

from __future__ import annotations

from .base import SurfaceScattering

__all__ = ["SMARTModel"]


class SMARTModel(SurfaceScattering):
    """Surface scattering placeholder for the SMART model."""

    MODEL_NAME = "SMART"
    __slots__ = ()

    def _compute_kirchhoff(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("SMARTModel is temporarily unavailable.")

    def _compute_complementary(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("SMARTModel is temporarily unavailable.")
