"""Example: Isotropic bare soil backscatter."""

import numpy as np

from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState

def main() -> None:
    config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)
    print("AIEM Backscatter Results:")
    for pol in (PolarizationState.HH, PolarizationState.VV, PolarizationState.HV):
        sigma = mwRTMs.compute_soil_backscatter(
            model="aiem",
            radar_config=config,
            frequency_ghz=5.4,
            rms_height_cm=1.0,
            correlation_length_cm=5.0,
            soil_permittivity=12.0 + 3.0j,
            polarization=pol,
        )
        print(f"  σ⁰_{pol.value.upper()} = {10.0 * np.log10(sigma + 1e-30):.2f} dB")


if __name__ == "__main__":
    main()
