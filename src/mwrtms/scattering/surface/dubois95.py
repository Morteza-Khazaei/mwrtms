"""Dubois95 empirical soil backscatter wrapper."""

from __future__ import annotations

from typing import Dict

from ssrt.surface.dubois95 import Dubois95

from ...medium.base import Medium
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["Dubois95Model"]


class Dubois95Model(PySSRTSurfaceScattering):
    """Expose the Dubois95 empirical model through the mwRTMs interface."""

    MODEL_NAME = "Dubois95"
    __slots__ = ()

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        model = Dubois95(
            fGHz=self.wave.frequency_ghz,
            theta=self.geometry.theta_i_deg,
            eps=permittivity,
            s=self.roughness.rms_height,
        )
        vv, hh, hv, vh = model.calc_sigma(todB=False)
        return {"vv": float(vv), "hh": float(hh), "hv": float(hv), "vh": float(vh)}
