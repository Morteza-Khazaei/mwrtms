"""Example: Vegetation backscatter coupling soil and canopy models."""

import numpy as np

from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState
from mwrtms.core import ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import SoilMedium, VegetationMedium
from mwrtms.scattering.volume import SSRTModel, CanopyProperties


def main() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    soil_medium = SoilMedium(
        moisture_m3m3=0.25,
        clay_fraction=0.3,
        sand_fraction=0.5,
    )
    perm = soil_medium.permittivity(frequency_hz=5.4e9)

    soil_sigma0 = {}
    for pol in (PolarizationState.HH, PolarizationState.VV, PolarizationState.HV):
        sigma = mwRTMs.compute_soil_backscatter(
            model="aiem",
            radar_config=config,
            frequency_ghz=5.4,
            rms_height_cm=1.0,
            correlation_length_cm=5.0,
            soil_permittivity=perm,
            polarization=pol,
        )
        soil_sigma0[pol.value] = sigma

    wave = ElectromagneticWave(frequency_hz=5.4e9)
    geometry = ScatteringGeometry(theta_i_deg=40.0)

    canopy = CanopyProperties(lai=3.0, vegetation_water_content=2.5, height=1.5)
    ssrt = SSRTModel(wave, geometry, canopy, soil_sigma0)

    veg_medium = VegetationMedium(gravimetric_moisture=0.6)
    air_medium = SoilMedium(moisture_m3m3=0.01, clay_fraction=0.3, sand_fraction=0.5)

    sigma_vv = ssrt.compute(air_medium, veg_medium, PolarizationState.VV)

    print("Vegetation (SSRT + AIEM):")
    print(f"  σ⁰_VV = {10 * np.log10(sigma_vv + 1e-30):.2f} dB")


if __name__ == "__main__":
    main()
