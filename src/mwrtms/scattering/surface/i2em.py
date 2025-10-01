"""Integral Equation Model 2 (I2EM) wrapper using pySSRT."""

from __future__ import annotations

from typing import Dict

from ssrt.surface.i2em import I2EM_Bistat_model

from ...medium.base import Medium
from ...utils.pyssrt import db_to_power
from ._pyssrt_bridge import PySSRTSurfaceScattering

__all__ = ["I2EMModel"]


class I2EMModel(PySSRTSurfaceScattering):
    """Expose the bistatic I2EM surface solver via the mwRTMs interface."""

    MODEL_NAME = "I2EM"
    __slots__ = ()

    def _run_pyssrt_model(self, medium_above: Medium, medium_below: Medium) -> Dict[str, float]:
        acf_label, power = self._acf_descriptor()
        sp_map = {"exp": 1, "gauss": 2, "pow": 3}
        sp = sp_map.get(acf_label, 1)
        xx = float(power) if (power is not None and sp == 3) else (1.5 if sp == 3 else 0.0)

        permittivity = medium_below.permittivity(self.wave.frequency_hz)
        theta_i = self.geometry.theta_i_deg
        theta_s = self.geometry.theta_s_deg
        phi_rel = (self.geometry.phi_s_deg - self.geometry.phi_i_deg) % 360.0

        vv_db, hh_db, hv_db, vh_db = I2EM_Bistat_model(
            fr=self.wave.frequency_ghz,
            sig=self.roughness.rms_height,
            L=self.roughness.correlation_length,
            thi=theta_i,
            ths=theta_s,
            phs=phi_rel,
            er=float(permittivity.real),
            sp=sp,
            xx=xx,
        )

        return {
            "vv": float(db_to_power(vv_db)),
            "hh": float(db_to_power(hh_db)),
            "hv": float(db_to_power(hv_db)),
            "vh": float(db_to_power(vh_db)),
        }
