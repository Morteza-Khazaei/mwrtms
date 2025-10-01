"""Advanced Integral Equation Method (AIEM) wrapper using pySSRT implementation."""

from __future__ import annotations

from typing import Dict

from ...pyssrt.surface.aiem import AIEM

from ...medium.base import Medium
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["AIEMModel"]


class AIEMModel(PySSRTSurfaceScattering):
    """Expose the legacy pySSRT AIEM solver through the mwRTMs interface."""

    MODEL_NAME = "AIEM"
    __slots__ = ("_include_multiple",)

    def __init__(self, wave, geometry, surface_roughness, *, include_multiple: bool = True, acf_override: str | None = None) -> None:
        super().__init__(wave, geometry, surface_roughness, acf_override=acf_override)
        self._include_multiple = bool(include_multiple)

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        acf_label, _ = self._acf_descriptor()
        surface_map = {"gauss": 1, "exp": 2, "pow": 3}
        surface_type = surface_map.get(acf_label, 2)

        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        theta_i = self.geometry.theta_i_deg
        theta_s = self.geometry.theta_s_deg
        phi_rel = (self.geometry.phi_s_deg - self.geometry.phi_i_deg) % 360.0

        k = self.wave.wavenumber
        kl = float(k * self.roughness.correlation_length)
        ks = float(k * self.roughness.rms_height)

        hh, vv, hv, vh = AIEM(
            theta_i=theta_i,
            theta_s=theta_s,
            phi_s=phi_rel,
            k=k,
            kl=kl,
            ks=ks,
            err=float(permittivity.real),
            eri=float(permittivity.imag),
            itype=surface_type,
            addMultiple=self._include_multiple,
            output_unit="linear",
        )

        return {"hh": float(hh), "vv": float(vv), "hv": float(hv), "vh": float(vh)}
