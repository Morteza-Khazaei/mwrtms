"""Bridges between mwRTMs abstractions and legacy pySSRT surface solvers."""

from __future__ import annotations

from typing import Dict, Tuple

from ...core.polarization import PolarizationState
from ...medium.base import Medium
from ...utils.pyssrt import determine_acf_descriptor, ensure_polarization_dict
from .base import SurfaceScattering

__all__ = ["PySSRTSurfaceScattering"]


class PySSRTSurfaceScattering(SurfaceScattering):
    """Base class that wraps a pySSRT surface backscatter implementation."""

    __slots__ = ("_acf_override", "_last_signature", "_sigma_cache")

    def __init__(self, wave, geometry, surface_roughness, *, acf_override: str | None = None) -> None:
        super().__init__(wave, geometry, surface_roughness)
        self._acf_override = acf_override
        self._last_signature: Tuple[float, ...] | None = None
        self._sigma_cache: Dict[str, float] | None = None

    def _compute_scattering(
        self, medium_above: Medium, medium_below: Medium, polarization: PolarizationState
    ) -> float:
        sigma = self._evaluate_sigma0(medium_above, medium_below)
        return float(sigma.get(polarization.value.lower(), 0.0))

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _evaluate_sigma0(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        signature = self._signature(medium_above, medium_below)
        if self._sigma_cache is None or signature != self._last_signature:
            sigma = self._run_pyssrt_model(medium_above, medium_below)
            self._sigma_cache = ensure_polarization_dict(sigma)
            self._last_signature = signature
        return self._sigma_cache

    def _signature(self, medium_above: Medium, medium_below: Medium) -> Tuple[float, ...]:
        eps_above = medium_above.permittivity(self.wave.frequency_hz)
        eps_below = medium_below.permittivity(self.wave.frequency_hz)
        return (
            float(self.wave.frequency_hz),
            float(self.geometry.theta_i_deg),
            float(self.geometry.theta_s_deg),
            float(self.geometry.phi_i_deg),
            float(self.geometry.phi_s_deg),
            float(eps_above.real),
            float(eps_above.imag),
            float(eps_below.real),
            float(eps_below.imag),
            float(self.roughness.rms_height),
            float(self.roughness.correlation_length),
        )

    def _acf_descriptor(self) -> Tuple[str, float | None]:
        return determine_acf_descriptor(self.roughness, self._acf_override)

    # ------------------------------------------------------------------
    # Hooks for subclasses
    # ------------------------------------------------------------------
    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        raise NotImplementedError
