import numpy as np

from mwrtms.scattering.surface.iem.guardrails import (
    diagnose_vka_complementary_slope_alignment,
)
from mwrtms.scattering.surface.iem.kirchhoff import VKA


def test_vka_complementary_slope_diagnostic_returns_results() -> None:
    theta_i = np.deg2rad(30.0)
    theta_s = np.deg2rad(20.0)
    phi_i = 0.0
    phi_s = 0.6

    diagnostics = diagnose_vka_complementary_slope_alignment(
        theta_i,
        theta_s,
        phi_i,
        phi_s,
        emit_warning=False,
    )

    assert "vka" in diagnostics
    assert "branches" in diagnostics

    branch_info = diagnostics["branches"]
    assert any(info["singular"] == 0.0 for info in branch_info.values())

    for info in branch_info.values():
        if info["singular"] == 0.0:
            assert np.isfinite(info["delta_norm"])


def test_vka_stationary_phase_slope_matches_tangent_backscatter() -> None:
    theta = np.deg2rad(30.0)
    vka = VKA(theta, theta, 0.0, np.pi, Rv=0.0 + 0.0j, Rh=0.0 + 0.0j)

    expected_slope = np.tan(theta)
    assert np.isclose(vka.zx, expected_slope, rtol=1e-9, atol=1e-9)
    assert np.isclose(vka.zy, 0.0, rtol=1e-9, atol=1e-9)
