"""Regression test comparing SPM3D predictions to the NMM3D LUT."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mwrtms.core import ElectromagneticWave, PolarizationState, ScatteringGeometry
from mwrtms.interface import ExponentialCorrelation, SurfaceRoughness
from mwrtms.medium import IsotropicMedium
from mwrtms.scattering.surface import SPM3DModel

_DATA_PATH = Path(__file__).resolve().parents[2] / "data" / "NMM3D_LUT_NRCS_40degree.dat"
_FREQUENCY_GHZ = 5.405
_INCIDENT_ANGLE_DEG = 40.0
_TEMPERATURE_K = 290.0


def _to_db(value: float) -> float:
    tiny = float(np.finfo(float).tiny)
    return float(10.0 * np.log10(max(value, tiny)))


def test_spm3d_matches_nmm3d_within_tolerance() -> None:
    table = np.loadtxt(_DATA_PATH)

    wave = ElectromagneticWave(_FREQUENCY_GHZ * 1e9)
    geometry = ScatteringGeometry(theta_i_deg=_INCIDENT_ANGLE_DEG)
    air = IsotropicMedium(1.0 + 0.0j, _TEMPERATURE_K)

    lambda_m = wave.wavelength

    vv_errors: list[float] = []
    hh_errors: list[float] = []
    hv_errors: list[float] = []

    for row in table:
        ratio = float(row[1])
        eps_r = float(row[2])
        eps_i = float(row[3])
        rms_norm = float(row[4])
        vv_ref = float(row[5])
        hh_ref = float(row[6])
        hv_ref = float(row[7])

        sigma = rms_norm * lambda_m
        corr_length = ratio * sigma
        roughness = SurfaceRoughness(
            rms_height=sigma,
            correlation_length=corr_length,
            correlation_function=ExponentialCorrelation(corr_length),
        )

        soil = IsotropicMedium(complex(eps_r, eps_i), _TEMPERATURE_K)
        model = SPM3DModel(wave, geometry, roughness)

        sigma_vv = model.compute(air, soil, PolarizationState.VV)
        sigma_hh = model.compute(air, soil, PolarizationState.HH)
        sigma_hv = model.compute(air, soil, PolarizationState.HV)

        vv_errors.append(_to_db(sigma_vv) - vv_ref)
        hh_errors.append(_to_db(sigma_hh) - hh_ref)
        if np.isfinite(hv_ref):
            hv_errors.append(_to_db(sigma_hv) - hv_ref)

    vv_errors_arr = np.asarray(vv_errors)
    hh_errors_arr = np.asarray(hh_errors)
    hv_errors_arr = np.asarray(hv_errors)

    vv_rmse = float(np.sqrt(np.mean(vv_errors_arr**2)))
    hh_rmse = float(np.sqrt(np.mean(hh_errors_arr**2)))
    hv_rmse = float(np.sqrt(np.mean(hv_errors_arr**2))) if hv_errors_arr.size else 0.0

    assert vv_rmse < 2.2, f"VV RMSE too large: {vv_rmse:.2f} dB"
    assert hh_rmse < 1.5, f"HH RMSE too large: {hh_rmse:.2f} dB"
    if hv_errors_arr.size:
        assert hv_rmse < 4.0, f"HV RMSE too large: {hv_rmse:.2f} dB"

    assert np.max(np.abs(vv_errors_arr)) < 5.0
    assert np.max(np.abs(hh_errors_arr)) < 4.0
