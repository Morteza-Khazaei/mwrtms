"""Example: Vegetation backscatter coupling soil and canopy models."""

import numpy as np

from mwrtms import mwRTMs
from mwrtms.core import ElectromagneticWave, ScatteringGeometry, PolarizationState
from mwrtms.medium import SoilMedium, VegetationMedium
from mwrtms.interface import ExponentialCorrelation, IsotropicRoughness
from mwrtms.scattering.volume import SSRTModel, CanopyProperties


def main() -> None:
    soil_result = mwRTMs.compute_soil_backscatter(
        "aiem",
        5.4,
        40.0,
        1.0,
        5.0,
        "exponential",
        soil_moisture=0.25,
    )

    soil_sigma0 = {pol: 10 ** (val / 10.0) for pol, val in soil_result.items()}

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
