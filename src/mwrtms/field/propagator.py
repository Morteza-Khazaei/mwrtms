"""Placeholder field propagator definitions."""

from __future__ import annotations

from dataclasses import dataclass

from ..core.geometry import ScatteringGeometry

__all__ = ["FieldPropagator"]


@dataclass
class FieldPropagator:
    """Simple container for propagator parameters."""

    u: float
    v: float
    q: float
    geometry: ScatteringGeometry

    def phase(self) -> float:
        return self.u * self.geometry.theta_i + self.v * self.geometry.theta_s
