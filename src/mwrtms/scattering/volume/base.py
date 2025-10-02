"""Base classes for volume scattering mechanisms."""

from __future__ import annotations

from ..base import ScatteringMechanism

__all__ = ["VolumeScattering"]


class VolumeScattering(ScatteringMechanism):
    """Abstract base for vegetation, snow, or atmospheric scattering."""

    pass
