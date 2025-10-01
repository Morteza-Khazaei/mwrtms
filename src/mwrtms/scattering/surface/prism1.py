"""PRISM-1 rough-surface backscatter wrapper."""

from __future__ import annotations

from typing import Dict

from ...pyssrt.surface.prism1 import PRISM1

from ...medium.base import Medium
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["PRISM1Model"]


class PRISM1Model(PySSRTSurfaceScattering):
    """Expose the pySSRT PRISM1 solver within the mwRTMs abstractions."""

    MODEL_NAME = "PRISM1"
    __slots__ = ()

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        model = PRISM1(
            f=self.wave.frequency_ghz,
            theta_i=self.geometry.theta_i_deg,
            eps=permittivity,
            s=self.roughness.rms_height,
        )
        vv, hh, hv, vh = model.calc_sigma(todB=False)
        return {"vv": float(vv), "hh": float(hh), "hv": float(hv), "vh": float(vh)}
