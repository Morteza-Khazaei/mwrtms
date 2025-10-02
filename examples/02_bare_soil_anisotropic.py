"""Example: Anisotropic bare soil (tilled rows)."""

import numpy as np

from mwrtms import ElectromagneticWave, RadarConfigurationFactory, PolarizationState
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.factory import ScatteringModelFactory
from mwrtms.medium import Medium


class _ConstantMedium(Medium):
    def __init__(self, permittivity: complex) -> None:
        super().__init__(temperature_k=293.15)
        self._eps = permittivity

    def permittivity(self, frequency_hz: float) -> complex:
        return self._eps


def main() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    wave = ElectromagneticWave(frequency_hz=5.4e9)
    surface = build_surface_from_statistics(
        rms_height_m=0.015,
        correlation_length_m=0.15,
        correlation_length_y_m=0.03,
    )

    model = ScatteringModelFactory.create_with_radar_config(
        "spm",
        config=config,
        wave=wave,
        surface=surface,
    )

    air = _ConstantMedium(1.0 + 0.0j)
    soil = _ConstantMedium(9.0 + 1.8j)

    sigma_vv = model.compute_with_config(air, soil, PolarizationState.VV, config)
    sigma_hh = model.compute_with_config(air, soil, PolarizationState.HH, config)

    print("Anisotropic (tilled soil):")
    print(f"  σ⁰_VV = {10.0 * np.log10(sigma_vv + 1e-30):.2f} dB")
    print(f"  σ⁰_HH = {10.0 * np.log10(sigma_hh + 1e-30):.2f} dB")
    print(f"  Anisotropy ratio ≈ {0.15 / 0.03:.1f}")


if __name__ == "__main__":
    main()
