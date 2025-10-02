"""Protocols for surface correlation functions."""

from __future__ import annotations

from typing import Protocol

__all__ = ["CorrelationFunction"]


class CorrelationFunction(Protocol):
    """Interface for roughness correlation functions.

    Implementations supply the roughness spectrum :math:`W^{(n)}` required by
    the surface scattering solvers. Using a ``Protocol`` enables duck typing so
    user-defined correlation models can be plugged in without altering the
    class hierarchy.
    """

    def spectrum(self, n: int, kx: float, ky: float, correlation_length: float) -> float:
        """Return the roughness spectrum ``W^{(n)}(k_x, k_y)``."""

        ...
