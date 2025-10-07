"""Utilities for IEM-family scattering models."""

from .correlation import (
    single_scale_spectrum_weights,
    multiscale_gaussian_weights,
    rms_slope_from_correlation,
    auto_spectral_order,
)
from .base import IEMBase, SurfaceRoughnessParameters
from .fresnel_utils import (
    compute_fresnel_incident,
    compute_fresnel_specular,
    compute_fresnel_nadir,
)
from .geometry_utils import (
    compute_q_vectors,
    compute_slope_components,
    compute_spatial_frequency,
)
from .spectrum_aiem import compute_aiem_spectrum
from .transition import compute_transition_function
from .kirchhoff import compute_kirchhoff_coefficients
from .complementary import (
    compute_expal,
    compute_complementary_vv,
    compute_complementary_hh,
    compute_complementary_hv,
    compute_complementary_vh,
)
from .i2em import I2EMModel
from .aiem import AIEMModel

__all__ = [
    "single_scale_spectrum_weights",
    "multiscale_gaussian_weights",
    "rms_slope_from_correlation",
    "auto_spectral_order",
    "IEMBase",
    "SurfaceRoughnessParameters",
    "compute_fresnel_incident",
    "compute_fresnel_specular",
    "compute_fresnel_nadir",
    "compute_q_vectors",
    "compute_slope_components",
    "compute_spatial_frequency",
    "compute_aiem_spectrum",
    "compute_transition_function",
    "compute_kirchhoff_coefficients",
    "compute_expal",
    "compute_complementary_vv",
    "compute_complementary_hh",
    "compute_complementary_hv",
    "compute_complementary_vh",
    "I2EMModel",
    "AIEMModel",
]
