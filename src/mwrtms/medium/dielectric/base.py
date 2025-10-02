"""Abstract base class for dielectric models."""

from __future__ import annotations

from abc import ABC, abstractmethod

__all__ = ["DielectricModel"]


class DielectricModel(ABC):
    """Common interface for dielectric calculations."""

    @abstractmethod
    def compute(self, frequency_hz: float, /, **kwargs) -> complex:
        """Return the complex relative permittivity ε_r(f)."""

        raise NotImplementedError
