"""Small Perturbation Method (SPM) surface model."""

from __future__ import annotations

import numpy as np

from ...factory import register_model
from .base import SurfaceScattering

__all__ = ["SPMModel"]


@register_model("spm")
class SPMModel(SurfaceScattering):
    """First-order SPM valid for kσ < 0.3."""

    MODEL_NAME = "SPM"

    def compute(self, medium_above, medium_below, polarization) -> float:
        R_h, R_v = self._compute_fresnel(medium_above, medium_below)
        return self._spm_backscatter(medium_above, medium_below, polarization, R_h, R_v)

    def _compute_kirchhoff(self, R_h, R_v, polarization):  # pragma: no cover - unused
        return 0.0

    def _compute_complementary(self, R_h, R_v, polarization):  # pragma: no cover - unused
        return 0.0

    def _spm_backscatter(self, medium_above, medium_below, polarization, R_h, R_v) -> float:
        theta = self._geometry.theta_i_rad
        sigma = self._surface.rms_height()
        ell = self._surface.correlation_length()

        eps1 = medium_above.permittivity(self._wave.frequency_hz)
        eps2 = medium_below.permittivity(self._wave.frequency_hz)

        k0 = self._wave.wavenumber
        k_inc = k0 * np.sqrt(eps1)
        k_tran = k0 * np.sqrt(eps2)

        sin_t = np.sin(theta)
        cos_t = np.cos(theta)

        kzi = k_inc * cos_t
        kxi = k_inc * sin_t
        k1zi = np.sqrt(k_tran**2 - kxi**2)

        kz = k_inc * cos_t
        kx = -k_inc * sin_t
        k1z = np.sqrt(k_tran**2 - kx**2)

        term_h = np.abs((kzi - k1zi) / (kzi + k1zi)) ** 2
        numerator_v = (k_tran**2 - k_inc**2) * (k_tran**2 * k_inc**2 * sin_t**2 + k_inc**2 * k1z * k1z)
        denominator_v = (k_tran**2 * kz + k_inc**2 * k1z) ** 2
        term_v = np.abs(numerator_v / denominator_v) ** 2

        prefactor = (
            8.0
            * k_inc**4
            * sigma**2
            * ell**2
            / (4.0 * k_inc**2 * ell**2 * sin_t**2 + 1.0) ** 1.5
            * cos_t**4
        )

        sigma_vv = np.real(prefactor * term_v)
        sigma_hh = np.real(prefactor * term_h)

        sqrt_eps = np.sqrt(max(np.real(eps2), 0.0))
        if sqrt_eps + 1.0 > 0.0:
            ratio_term = max((sqrt_eps - 1.0) / (sqrt_eps + 1.0), 0.0)
        else:
            ratio_term = 0.0
        hv_coeff = 0.23 * np.sqrt(ratio_term)
        krms = k_inc * sigma
        sigma_hv = hv_coeff * (1.0 - np.exp(-krms)) * sigma_vv

        if polarization.value in ("vv", "v"):
            return float(max(sigma_vv, 0.0))
        if polarization.value in ("hh", "h"):
            return float(max(sigma_hh, 0.0))
        if polarization.value in ("hv", "vh"):
            return float(max(sigma_hv, 0.0))
        return 0.0
