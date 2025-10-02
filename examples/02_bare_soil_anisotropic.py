"""Example: Anisotropic bare soil (tilled rows)."""

import numpy as np

from mwrtms import ElectromagneticWave, PolarizationState, RadarConfigurationFactory, ScatteringScene
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.medium import HomogeneousMedium


def main() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    wave = ElectromagneticWave(frequency_hz=5.4e9)
    surface = build_surface_from_statistics(
        rms_height_m=0.015,
        correlation_length_m=0.15,
        correlation_length_y_m=0.03,
    )

    scene = ScatteringScene(
        radar_config=config,
        wave=wave,
        medium_above=HomogeneousMedium(1.0 + 0.0j),
        medium_below=HomogeneousMedium(9.0 + 1.8j),
        surface=surface,
    )

    result = scene.run_model("spm")

    print("Anisotropic (tilled soil):")
    print(f"  σ⁰_VV = {10.0 * np.log10(result[PolarizationState.VV] + 1e-30):.2f} dB")
    print(f"  σ⁰_HH = {10.0 * np.log10(result[PolarizationState.HH] + 1e-30):.2f} dB")
    print(f"  Anisotropy ratio ≈ {0.15 / 0.03:.1f}")


if __name__ == "__main__":
    main()
