"""Small Perturbation Method (SPM3D) wrapper."""

from __future__ import annotations

from typing import Dict

from ssrt.surface.spm import SPM3D

from ...medium.base import Medium
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["SPM3DModel"]


class SPM3DModel(PySSRTSurfaceScattering):
    """Expose the pySSRT SPM3D solver in the mwRTMs interface."""

    MODEL_NAME = "SPM3D"
    __slots__ = ()

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        model = SPM3D(
            fr=self.wave.frequency_ghz,
            sig=self.roughness.rms_height,
            L=self.roughness.correlation_length,
            thi=self.geometry.theta_i_deg,
            eps=permittivity,
        )
        vv, hh, hv, vh = model.calc_sigma(todB=False)
        return {"vv": float(vv), "hh": float(hh), "hv": float(hv), "vh": float(vh)}
