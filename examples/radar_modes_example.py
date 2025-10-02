"""Demonstrate radar configuration helpers and scattering workflows."""

from __future__ import annotations

import numpy as np

from mwrtms import (
    ElectromagneticWave,
    PolarizationState,
    RadarConfigurationFactory,
    RadarConfiguration,
    ScatteringGeometry,
    ScatteringScene,
    mwRTMsFacade,
)
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics


def _print_result(label: str, sigma: float) -> None:
    print(f"{label}: {10.0 * np.log10(sigma + 1e-30):.2f} dB")


def monostatic_example() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    wave = ElectromagneticWave(frequency_hz=5.405e9)
    surface = build_surface_from_statistics(0.008, 0.04)
    scene = ScatteringScene(
        radar_config=config,
        wave=wave,
        medium_above=HomogeneousMedium(1.0 + 0.0j),
        medium_below=HomogeneousMedium(12.0 + 2.5j),
        surface=surface,
    )
    result = scene.run_model("spm")
    _print_result("Monostatic σ₀(VV)", result[PolarizationState.VV])


def bistatic_example() -> None:
    config = RadarConfigurationFactory.create_bistatic(
        theta_i_deg=30.0,
        theta_s_deg=45.0,
        phi_s_deg=90.0,
    )
    wave = ElectromagneticWave(frequency_hz=5.405e9)
    surface = build_surface_from_statistics(0.005, 0.03)

    scene = ScatteringScene(
        radar_config=config,
        wave=wave,
        medium_above=HomogeneousMedium(1.0 + 0.0j),
        medium_below=HomogeneousMedium(8.0 + 1.5j),
        surface=surface,
    )
    result = scene.run_model("spm")
    sigma_vv = result[PolarizationState.VV]
    angle = config.bistatic_angle()
    _print_result(f"Bistatic σ₀(VV) (angle {angle:.1f}°)", sigma_vv)


def multi_angle_example() -> None:
    configs = RadarConfigurationFactory.create_multi_angle_monostatic([20.0, 30.0, 40.0])
    wave = ElectromagneticWave(frequency_hz=1.0e9)
    surface = build_surface_from_statistics(0.004, 0.02)

    for cfg in configs:
        scene = ScatteringScene(
            radar_config=cfg,
            wave=wave,
            medium_above=HomogeneousMedium(1.0 + 0.0j),
            medium_below=HomogeneousMedium(5.0 + 0.5j),
            surface=surface,
        )
        sigma_hh = scene.run_scalar("spm", PolarizationState.HH)
        _print_result(f"Multi-angle θ={cfg.look_angle_deg:.1f}° σ₀(HH)", sigma_hh)


def manual_vs_factory() -> None:
    # Manual geometry construction
    manual_geometry = ScatteringGeometry(theta_i_deg=35.0, theta_s_deg=35.0)
    manual_config = RadarConfigurationFactory.from_geometry(manual_geometry, description="Manual geometry")

    factory_config = RadarConfigurationFactory.create_monostatic(theta_deg=35.0, description="Factory geometry")
    print(
        "Manual vs factory geometry match:",
        manual_config.matches_geometry(factory_config.geometry),
    )


def facade_example() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    sigma_vv = mwRTMsFacade.compute_soil_backscatter(
        model="spm",
        radar_config=config,
        frequency_ghz=5.4,
        rms_height_cm=0.8,
        correlation_length_cm=4.0,
        soil_permittivity=10.0 + 2.0j,
        correlation="exponential",
        polarization=PolarizationState.VV,
    )
    _print_result("Facade σ₀(VV)", sigma_vv)


def compute_with_config_usage() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=25.0)
    wave = ElectromagneticWave(frequency_hz=3.0e9)
    surface = build_surface_from_statistics(0.002, 0.015)

    scene = ScatteringScene(
        radar_config=config,
        wave=wave,
        medium_above=HomogeneousMedium(1.0 + 0.0j),
        medium_below=HomogeneousMedium(6.0 + 0.9j),
        surface=surface,
    )
    sigma_hv = scene.run_scalar("spm", PolarizationState.HV)
    _print_result("Compute-with-config σ₀(HV)", sigma_hv)


if __name__ == "__main__":
    print("--- Monostatic example ---")
    monostatic_example()

    print("\n--- Bistatic example ---")
    bistatic_example()

    print("\n--- Multi-angle monostatic example ---")
    multi_angle_example()

    print("\n--- Manual geometry vs factory ---")
    manual_vs_factory()

    print("\n--- Facade helper ---")
    facade_example()

    print("\n--- Direct compute_with_config usage ---")
    compute_with_config_usage()
