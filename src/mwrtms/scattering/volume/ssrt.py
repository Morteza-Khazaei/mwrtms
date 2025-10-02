"""Placeholder SSRT volume model while implementation is unavailable."""

from __future__ import annotations

from ...medium.base import Medium
from .base import VolumeScattering

__all__ = ["SSRTModel"]


class SSRTModel(VolumeScattering):
    """Volume scattering placeholder for the SSRT model."""

    MODEL_NAME = "SSRT"
    __slots__ = ()

    def _compute_volume_response(self, medium_above: Medium, medium_below: Medium, polarization):  # pragma: no cover - placeholder
        raise NotImplementedError("SSRTModel is temporarily unavailable.")
