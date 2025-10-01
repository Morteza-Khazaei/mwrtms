"""Volume scattering base classes."""

from __future__ import annotations

from abc import ABC, abstractmethod

from ...core.polarization import PolarizationState
from ...medium.base import Medium
from ..base import ScatteringMechanism

__all__ = ["VolumeScattering"]


class VolumeScattering(ScatteringMechanism, ABC):
    """Base class for vegetation/atmospheric scattering models."""

    __slots__ = ()

    def _compute_scattering(self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState) -> float:
        return self._compute_volume_response(medium_above, medium_below, polarization)

    @abstractmethod
    def _compute_volume_response(
        self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState
    ) -> float:
        """Volume scattering contribution implemented by subclasses."""
