"""Composite surface generator."""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from ..base import Surface
from .base import SurfaceGenerator


class CompositeSurfaceGenerator(SurfaceGenerator):
    """Combine multiple surfaces into a weighted composite."""

    def __init__(
        self,
        surfaces: Sequence[Surface],
        weights: Sequence[float] | None = None,
    ) -> None:
        if not surfaces:
            raise ValueError("at least one surface required")

        shapes = {s.shape for s in surfaces}
        if len(shapes) != 1:
            raise ValueError("all surfaces must share the same shape")

        sizes = {s.physical_size for s in surfaces}
        if len(sizes) != 1:
            raise ValueError("all surfaces must share the same physical size")

        shape = next(iter(shapes))
        size = next(iter(sizes))
        super().__init__(shape, size)

        self._surfaces = list(surfaces)
        if weights is None:
            self._weights = [1.0] * len(self._surfaces)
        else:
            if len(weights) != len(self._surfaces):
                raise ValueError("weights must match number of surfaces")
            self._weights = list(float(w) for w in weights)

    def _generate_impl(self) -> Surface:
        heights = np.zeros(self._shape, dtype=float)
        for surface, weight in zip(self._surfaces, self._weights):
            heights += weight * np.asarray(surface)

        metadata = {
            "generator": "composite",
            "n_surfaces": len(self._surfaces),
            "weights": list(self._weights),
        }
        return Surface(heights, self._physical_size, metadata)
