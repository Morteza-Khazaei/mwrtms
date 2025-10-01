"""SMART semi-empirical soil backscatter wrapper."""

from __future__ import annotations

from typing import Dict

from ...pyssrt.surface.smart import SMART

from ...medium.base import Medium
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["SMARTModel"]


class SMARTModel(PySSRTSurfaceScattering):
    """Expose the pySSRT SMART model through mwRTMs abstractions."""

    MODEL_NAME = "SMART"
    __slots__ = ()

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        model = SMART(
            fGHz=self.wave.frequency_ghz,
            theta_deg=self.geometry.theta_i_deg,
            s=self.roughness.rms_height,
            eps=permittivity,
        )
        vv, hh, hv, vh = model.calc_sigma(todB=False)
        return {"vv": float(vv), "hh": float(hh), "hv": float(hv), "vh": float(vh)}
