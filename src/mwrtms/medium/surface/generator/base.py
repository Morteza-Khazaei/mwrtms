"""Base class for surface generators."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from ..base import Surface


class SurfaceGenerator:
    """Template for surface generators with optional random seeding."""

    def __init__(self, shape: Tuple[int, int], physical_size: Tuple[float, float]) -> None:
        self._shape = shape
        self._physical_size = physical_size
        self._random_seed: Optional[int] = None

    def generate(self) -> Surface:
        """Generate a new :class:`Surface` instance."""

        self._initialize()
        surface = self._generate_impl()
        return self._post_process(surface)

    def _initialize(self) -> None:
        if self._random_seed is not None:
            np.random.seed(self._random_seed)

    def _post_process(self, surface: Surface) -> Surface:
        """Standard post-processing pipeline (centre the surface)."""

        return surface.center()

    def setRandomSeed(self, seed: int) -> "SurfaceGenerator":  # noqa: N802 - Tamaas-inspired API
        """Set the RNG seed (fluent API)."""

        self._random_seed = seed
        return self

    def _generate_impl(self) -> Surface:
        raise NotImplementedError
