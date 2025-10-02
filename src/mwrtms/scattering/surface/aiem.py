"""Placeholder AIEM model keeping interface while implementation is unavailable."""

from __future__ import annotations

from .base import SurfaceScattering

__all__ = ["AIEMModel"]


class AIEMModel(SurfaceScattering):
    """Surface scattering placeholder for the AIEM model."""

    MODEL_NAME = "AIEM"
    __slots__ = ()

    def _compute_kirchhoff(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("AIEMModel is temporarily unavailable.")

    def _compute_complementary(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("AIEMModel is temporarily unavailable.")
