"""Small Perturbation Method (SPM3D) surface backscatter model."""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np

from ...core.polarization import PolarizationState
from ...interface.roughness import SurfaceRoughness
from ...medium.base import Medium
from .base import SurfaceScattering

__all__ = ["SPM3DModel"]


class SPM3DModel(SurfaceScattering):
    """First-order small perturbation method (SPM3D) implementation.

    The class adheres to the mwRTMs object model by:

    * **Encapsulation** – the surface descriptors and cached results are private
      and exposed only via the public `compute` method inherited from the base
      class.
    * **Inheritance** – `SurfaceScattering` supplies the common interface;
      SPM3D overrides the scattering evaluation with the analytical solution.
    * **Polymorphism** – the model can be created through the factory/facade and
      used interchangeably with other surface mechanisms.
    * **Abstraction** – solver details are confined to private helpers, keeping
      the external API compact.
    """

    MODEL_NAME = "SPM3D"
    __slots__ = ("_include_cross_pol", "_sigma_cache", "_last_signature")
    _mutable_slots = {"_sigma_cache", "_last_signature"}

    def __init__(
        self,
        wave,
        geometry,
        surface_roughness: SurfaceRoughness,
        *,
        include_cross_pol: bool = True,
    ) -> None:
        super().__init__(wave, geometry, surface_roughness)
        self._include_cross_pol = bool(include_cross_pol)
        self._sigma_cache: Dict[str, float] | None = None
        self._last_signature: Tuple[float, ...] | None = None

    # ------------------------------------------------------------------
    # SurfaceScattering API
    # ------------------------------------------------------------------
    def _compute_scattering(
        self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState
    ) -> float:
        sigma = self._evaluate_sigma0(medium_above, medium_below)
        return float(sigma.get(polarization.value.lower(), 0.0))

    def _compute_kirchhoff(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        # The classical SPM3D solution is implemented via `_compute_scattering`.
        # The abstract method is still provided to satisfy the inheritance contract.
        return 0.0

    def _compute_complementary(self, fresnel: dict[str, complex], polarization: PolarizationState) -> float:
        return 0.0

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _evaluate_sigma0(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        signature = self._signature(medium_below)
        if self._sigma_cache is None or signature != self._last_signature:
            self._sigma_cache = self._compute_sigma0(medium_above, medium_below)
            self._last_signature = signature
        return self._sigma_cache

    def _signature(self, medium_below: Medium) -> Tuple[float, ...]:
        eps = medium_below.permittivity(self.wave.frequency_hz)
        return (
            float(self.wave.frequency_hz),
            float(self.geometry.theta_i_deg),
            float(self.roughness.rms_height),
            float(self.roughness.correlation_length),
            float(np.real(eps)),
            float(np.imag(eps)),
            float(self._include_cross_pol),
        )

    def _compute_sigma0(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        theta = self.geometry.theta_i
        h = self.roughness.rms_height
        corr_length = self.roughness.correlation_length
        k = self.wave.wavenumber
        eps = complex(medium_below.permittivity(self.wave.frequency_hz))

        k1 = k * np.sqrt(eps)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        kzi = k * cos_theta
        kxi = k * sin_theta
        k1zi = np.sqrt(k1**2 - kxi**2)

        kz = k * cos_theta
        kx = -k * sin_theta
        k1z = np.sqrt(k1**2 - kx**2)

        term_h = abs((kzi - k1zi) / (kzi + k1zi)) ** 2

        numerator_v = (k1**2 - k**2) * (k1**2 * k**2 * sin_theta**2 + k**2 * k1z * k1z)
        denominator_v = (k1**2 * kz + k**2 * k1z) ** 2
        term_v = abs(numerator_v / denominator_v) ** 2

        sin_theta_sq = sin_theta**2
        denom = np.power(4.0 * k**2 * corr_length**2 * sin_theta_sq + 1.0, 1.5)
        prefactor = 8.0 * (k**4) * (h**2) * (corr_length**2) * cos_theta**4 / denom

        sig_hh = float(np.real(prefactor * term_h))
        sig_vv = float(np.real(prefactor * term_v))

        tiny = float(np.finfo(float).tiny)
        sig_hh = max(sig_hh, tiny)
        sig_vv = max(sig_vv, tiny)

        if self._include_cross_pol:
            eps_real = max(float(np.real(eps)), 0.0)
            sqrt_eps_real = float(np.sqrt(eps_real)) if eps_real > 0 else 0.0
            if sqrt_eps_real > 0:
                ratio_term = (sqrt_eps_real - 1.0) / (sqrt_eps_real + 1.0)
                ratio_term = max(ratio_term, 0.0)
                hv_coeff = 0.23 * float(np.sqrt(ratio_term))
            else:
                hv_coeff = 0.0
            krms = k * h
            sig_hv = hv_coeff * (1.0 - float(np.exp(-krms))) * sig_vv
        else:
            sig_hv = 0.0

        sig_hv = max(float(sig_hv), tiny)

        return {
            "vv": sig_vv,
            "hh": sig_hh,
            "hv": sig_hv,
            "vh": sig_hv,
        }
