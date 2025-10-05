"""Debug script to examine actual HV values from AIEM vs NMM3D."""

import numpy as np
from pathlib import Path
from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState

# Load NMM3D data
lut_path = Path("data/NMM3D_LUT_NRCS_40degree.dat")
table = np.loadtxt(lut_path)

# Filter for 40 degrees
mask = np.isclose(table[:, 0], 40.0)
rows = table[mask]

# Take first few samples
print("=" * 80)
print("Comparing AIEM HV with NMM3D HV (first 10 samples)")
print("=" * 80)

frequency_ghz = 5.405
wavelength_m = 0.3 / frequency_ghz
radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)

for i, row in enumerate(rows[:10]):
    theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
    
    sigma = rms_norm * wavelength_m
    corr_len = ratio * sigma
    rms_height_cm = sigma * 100.0
    correlation_length_cm = corr_len * 100.0
    soil_permittivity = complex(float(eps_r), float(eps_i))
    
    # Compute AIEM HV
    hv_aiem = mwRTMs.compute_soil_backscatter(
        model='aiem',
        radar_config=radar_config,
        frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm,
        correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity,
        correlation='exponential',
        polarization=PolarizationState.HV,
    )
    
    hv_aiem_db = 10.0 * np.log10(hv_aiem) if hv_aiem > 0 else -999
    
    print(f"\nSample {i+1}:")
    print(f"  Params: ℓ/σ={ratio:.1f}, εr={eps_r:.1f}-j{eps_i:.1f}, kσ={2*np.pi*rms_norm:.3f}")
    print(f"  NMM3D HV: {hv_ref:7.2f} dB")
    print(f"  AIEM HV:  {hv_aiem_db:7.2f} dB (σ⁰={hv_aiem:.3e})")
    print(f"  Diff:     {hv_aiem_db - hv_ref:7.2f} dB")

print("\n" + "=" * 80)
