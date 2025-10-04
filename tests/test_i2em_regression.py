"""Regression test comparing I2EM predictions to the NMM3D LUT."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from mwrtms.core import ElectromagneticWave, PolarizationState, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface import I2EMModel

_DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "NMM3D_LUT_NRCS_40degree.dat"
_FREQUENCY_GHZ = 5.405
_INCIDENT_ANGLE_DEG = 40.0


def _to_db(value: float) -> float:
    """Safely convert a linear value to decibels."""
    tiny = float(np.finfo(float).tiny)
    return float(10.0 * np.log10(max(value, tiny)))


def test_i2em_matches_nmm3d_within_tolerance() -> None:
    """
    Validates the I2EM model against the NMM3D LUT.

    This test ensures that the co-polarization (VV, HH) results from the
    I2EMModel remain consistent with the reference NMM3D data, preventing
    future regressions. It calculates the Root Mean Square Error (RMSE)
    and asserts that it stays below a pre-defined threshold derived from
    successful validation runs.
    """
    table = np.loadtxt(_DATA_PATH)

    wave = ElectromagneticWave(_FREQUENCY_GHZ * 1e9)
    geometry = ScatteringGeometry(theta_i_deg=_INCIDENT_ANGLE_DEG)
    air = HomogeneousMedium(1.0 + 0.0j)
    lambda_m = wave.wavelength

    vv_errors: list[float] = []
    hh_errors: list[float] = []

    for row in table:
        _, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, _ = row

        sigma = rms_norm * lambda_m
        surface = build_surface_from_statistics(
            sigma,
            ratio * sigma,
            correlation_type="exponential",
        )
        soil = HomogeneousMedium(complex(eps_r, eps_i))
        model = I2EMModel(wave, geometry, surface)

        sigma_vv = model.compute(air, soil, PolarizationState.VV)
        sigma_hh = model.compute(air, soil, PolarizationState.HH)

        vv_errors.append(_to_db(sigma_vv) - vv_ref)
        hh_errors.append(_to_db(sigma_hh) - hh_ref)

    vv_errors_arr = np.asarray(vv_errors)
    hh_errors_arr = np.asarray(hh_errors)

    vv_rmse = float(np.sqrt(np.mean(vv_errors_arr**2)))
    hh_rmse = float(np.sqrt(np.mean(hh_errors_arr**2)))

    assert vv_rmse < 1.5, f"I2EM VV RMSE is too high: {vv_rmse:.2f} dB"
    assert hh_rmse < 1.0, f"I2EM HH RMSE is too high: {hh_rmse:.2f} dB"