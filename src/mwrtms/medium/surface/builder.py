"""Helpers to synthesize surfaces from analytic statistics."""

from __future__ import annotations

from functools import lru_cache
from typing import Optional

from .base import Surface
from .generator.synthetic import SyntheticSurfaceGenerator

_TEMPLATE_SHAPE = (256, 256)
_TEMPLATE_SIZE = (1.0, 1.0)


@lru_cache(maxsize=1)
def _template_surface(seed: int = 2024) -> Surface:
    generator = SyntheticSurfaceGenerator(_TEMPLATE_SHAPE, _TEMPLATE_SIZE).setRandomSeed(seed)
    surface = generator.generate()
    return surface.center()


def build_surface_from_statistics(
    rms_height_m: float,
    correlation_length_m: float,
    *,
    correlation_length_y_m: Optional[float] = None,
    correlation_type: str = "exponential",
) -> Surface:
    """Construct a synthetic :class:`Surface` matching basic statistics."""

    if rms_height_m <= 0.0:
        raise ValueError("rms_height_m must be positive")
    if correlation_length_m <= 0.0:
        raise ValueError("correlation_length_m must be positive")
    if correlation_length_y_m is not None and correlation_length_y_m <= 0.0:
        raise ValueError("correlation_length_y_m must be positive")

    template = _template_surface()

    surface = Surface(
        template.heights,
        template.physical_size,
        {**template.metadata, "correlation_type": correlation_type},
    )

    surface = surface.scale_rms(rms_height_m)

    if correlation_length_y_m is None:
        base_corr = max(surface.correlation_length(), 1e-12)
        scale = correlation_length_m / base_corr
        lx, ly = surface.physical_size
        surface = surface.with_physical_size((lx * scale, ly * scale))
    else:
        base_x, base_y = surface.correlation_lengths_xy()
        base_x = max(base_x, 1e-12)
        base_y = max(base_y, 1e-12)
        scale_x = correlation_length_m / base_x
        scale_y = correlation_length_y_m / base_y
        lx, ly = surface.physical_size
        surface = surface.with_physical_size((lx * scale_x, ly * scale_y))

    metadata = surface.metadata
    metadata["correlation_type"] = correlation_type
    metadata["correlation_length_target"] = correlation_length_m
    if correlation_length_y_m is not None:
        metadata["correlation_length_x_target"] = correlation_length_m
        metadata["correlation_length_y_target"] = correlation_length_y_m
    return Surface(surface.heights, surface.physical_size, metadata)
