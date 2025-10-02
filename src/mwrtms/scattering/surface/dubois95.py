"""Placeholder Dubois95 model while implementation is unavailable."""

from __future__ import annotations

from .base import SurfaceScattering

__all__ = ["Dubois95Model"]


class Dubois95Model(SurfaceScattering):
    """Surface scattering placeholder for the Dubois95 model."""

    MODEL_NAME = "Dubois95"
    __slots__ = ()

    def _compute_kirchhoff(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("Dubois95Model is temporarily unavailable.")

    def _compute_complementary(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("Dubois95Model is temporarily unavailable.")
