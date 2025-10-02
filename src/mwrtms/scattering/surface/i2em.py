"""Placeholder I2EM model while implementation is unavailable."""

from __future__ import annotations

from .base import SurfaceScattering

__all__ = ["I2EMModel"]


class I2EMModel(SurfaceScattering):
    """Surface scattering placeholder for the I2EM model."""

    MODEL_NAME = "I2EM"
    __slots__ = ()

    def _compute_kirchhoff(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("I2EMModel is temporarily unavailable.")

    def _compute_complementary(self, *args, **kwargs):  # pragma: no cover - placeholder
        raise NotImplementedError("I2EMModel is temporarily unavailable.")
