"""Statistical helpers for :class:`~mwrtms.medium.surface.base.Surface`."""

from __future__ import annotations

from typing import Dict

from .base import Surface


def summary(surface: Surface) -> Dict[str, float]:
    """Return a concise summary of common statistical parameters."""

    stats = {
        "mean_height": surface.mean_height(),
        "rms_height": surface.rms_height(),
        "peak_to_valley": surface.peak_to_valley(),
        "rms_slope_x": surface.rms_slope("x"),
        "rms_slope_y": surface.rms_slope("y"),
        "rms_slope_total": surface.rms_slope(),
    }

    if not surface.is_isotropic():
        ell_x, ell_y = surface.correlation_lengths_xy()
        stats["correlation_length_x"] = ell_x
        stats["correlation_length_y"] = ell_y
    else:
        stats["correlation_length"] = surface.correlation_length()

    return stats


def isotropy(surface: Surface, tolerance: float = 0.2) -> Dict[str, float]:
    """Return isotropy indicators for convenience."""

    return {
        "is_isotropic": surface.is_isotropic(tolerance=tolerance),
        "anisotropy_ratio": surface.anisotropy_ratio(),
    }
