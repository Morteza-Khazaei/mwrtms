"""Placeholder PRISM1 model while implementation is unavailable."""

from __future__ import annotations

from .base import SurfaceScattering

__all__ = ["PRISM1Model"]


class PRISM1Model(SurfaceScattering):
    """Surface scattering placeholder for the PRISM1 model."""

    MODEL_NAME = "PRISM1"
    __slots__ = ()

    def _compute_kirchhoff(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("PRISM1Model is temporarily unavailable.")

    def _compute_complementary(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("PRISM1Model is temporarily unavailable.")
