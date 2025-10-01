"""Single-scattering radiative transfer (SSRT) wrapper around pySSRT."""

from __future__ import annotations

from typing import Dict, Tuple

from ...pyssrt.core import S2RTR

from ...core.polarization import PolarizationState
from ...interface.roughness import SurfaceRoughness
from ...medium.base import Medium
from ...utils.pyssrt import determine_acf_descriptor, ensure_polarization_dict
from .base import VolumeScattering

__all__ = ["SSRTModel"]


class SSRTModel(VolumeScattering):
    """Couple soil and vegetation scattering using the legacy pySSRT solver."""

    MODEL_NAME = "SSRT"
    __slots__ = (
        "_roughness",
        "_single_scatter_albedo",
        "_extinction_coefficient",
        "_canopy_thickness",
        "_surface_model",
        "_canopy_model",
        "_acf_override",
        "_last_signature",
        "_sigma_cache",
    )

    _SURFACE_NORMALISATION = {
        "aiem": "AIEM",
        "prism1": "PRISM1",
        "dubois95": "Dubois95",
        "smart": "SMART",
        "spm3d": "SPM3D",
        "spm": "SPM3D",
        "i2em": "I2EM",
    }
    _CANOPY_NORMALISATION = {"diff": "Diff", "spec": "Spec"}

    def __init__(
        self,
        wave,
        geometry,
        surface_roughness: SurfaceRoughness,
        *,
        single_scatter_albedo: float,
        extinction_coefficient: float,
        canopy_thickness: float,
        surface_model: str = "I2EM",
        canopy_model: str = "Diff",
        acf_override: str | None = None,
    ) -> None:
        super().__init__(wave, geometry)
        if not 0.0 <= single_scatter_albedo < 1.0:
            raise ValueError("single_scatter_albedo must be within [0, 1)")
        if extinction_coefficient <= 0.0:
            raise ValueError("extinction_coefficient must be positive")
        if canopy_thickness <= 0.0:
            raise ValueError("canopy_thickness must be positive")

        try:
            surf_label = self._SURFACE_NORMALISATION[surface_model.lower()]
        except KeyError as exc:
            available = ", ".join(sorted(self._SURFACE_NORMALISATION.values()))
            raise ValueError(f"surface_model must be one of {available}") from exc

        try:
            canopy_label = self._CANOPY_NORMALISATION[canopy_model.lower()]
        except KeyError as exc:
            available = ", ".join(sorted(self._CANOPY_NORMALISATION.values()))
            raise ValueError(f"canopy_model must be one of {available}") from exc

        self._roughness = surface_roughness
        self._single_scatter_albedo = float(single_scatter_albedo)
        self._extinction_coefficient = float(extinction_coefficient)
        self._canopy_thickness = float(canopy_thickness)
        self._surface_model = surf_label
        self._canopy_model = canopy_label
        self._acf_override = acf_override
        self._last_signature: Tuple[float, ...] | None = None
        self._sigma_cache: Dict[str, float] | None = None

    # ------------------------------------------------------------------
    # VolumeScattering interface
    # ------------------------------------------------------------------
    def _compute_volume_response(
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
            acf_label, _ = determine_acf_descriptor(self._roughness, self._acf_override)
            ssrt = S2RTR(
                frq_GHz=self.wave.frequency_ghz,
                theta_i=self.geometry.theta_i_deg,
                theta_s=self.geometry.theta_s_deg,
                phi_i=self.geometry.phi_i_deg,
                phi_s=self.geometry.phi_s_deg,
                s=self._roughness.rms_height,
                cl=self._roughness.correlation_length,
                eps2=medium_above.permittivity(self.wave.frequency_hz),
                eps3=medium_below.permittivity(self.wave.frequency_hz),
                a=self._single_scatter_albedo,
                kappa_e=self._extinction_coefficient,
                d=self._canopy_thickness,
                acftype=acf_label,
                RT_models={"RT_s": self._surface_model, "RT_c": self._canopy_model},
            )
            result = ssrt.calc_sigma(todB=False)
            if isinstance(result, tuple):
                total = result[-1]
            else:
                total = result
            self._sigma_cache = ensure_polarization_dict(total)
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
            float(self._roughness.rms_height),
            float(self._roughness.correlation_length),
            float(eps_above.real),
            float(eps_above.imag),
            float(eps_below.real),
            float(eps_below.imag),
            self._single_scatter_albedo,
            self._extinction_coefficient,
            self._canopy_thickness,
        )
