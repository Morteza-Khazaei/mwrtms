"""Abstract base class for electromagnetic media."""

from __future__ import annotations

from abc import ABC, abstractmethod

__all__ = ["Medium"]


class Medium(ABC):
    """Abstract base class representing a homogeneous medium.

    The class showcases *abstraction* by declaring the interface shared by all
    media and *encapsulation* by keeping the temperature and dielectric cache
    private. Concrete subclasses (soil, vegetation, ocean, snow, …) provide the
    specific permittivity behaviour.
    """

    __slots__ = ("_temperature_k", "_dielectric_cache")

    def __init__(self, temperature_k: float = 293.15) -> None:
        if temperature_k <= 0.0:
            raise ValueError("temperature_k must be positive")
        self._temperature_k = float(temperature_k)
        self._dielectric_cache: dict[float, complex] = {}

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------
    @abstractmethod
    def permittivity(self, frequency_hz: float) -> complex:
        """Return the complex relative permittivity ε_r(f)."""

    # ------------------------------------------------------------------
    # Encapsulated accessors
    # ------------------------------------------------------------------
    @property
    def temperature_k(self) -> float:
        """Physical temperature of the medium in Kelvin (read-only)."""

        return self._temperature_k

    @property
    def permeability(self) -> complex:
        """Relative magnetic permeability μ_r (unity for natural media)."""

        return 1.0 + 0.0j
