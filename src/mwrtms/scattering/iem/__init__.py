"""Utilities for IEM-family scattering models."""

from .correlation import (
    single_scale_spectrum_weights,
    multiscale_gaussian_weights,
    rms_slope_from_correlation,
    auto_spectral_order,
)
from .base import IEMBase, SurfaceRoughnessParameters

__all__ = [
    "single_scale_spectrum_weights",
    "multiscale_gaussian_weights",
    "rms_slope_from_correlation",
    "auto_spectral_order",
    "IEMBase",
    "SurfaceRoughnessParameters",
]
