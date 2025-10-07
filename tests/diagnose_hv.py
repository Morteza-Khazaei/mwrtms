"""Diagnostic script to understand HV cross-pol issue."""

import numpy as np
from pathlib import Path
from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState

# Load NMM3D data
lut_path = Path("data/NMM3D_LUT_NRCS_40degree.dat")
table = np.loadtxt(lut_path)
mask = np.isclose(table[:, 0], 40.0)
rows = table[mask]

frequency_ghz = 5.405
wavelength_m = 0.3 / frequency_ghz
radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)

print("=" * 80)
print("HV Cross-Polarization Diagnostic")
print("=" * 80)

# Take first sample
row = rows[0]
theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row

sigma = rms_norm * wavelength_m
corr_len = ratio * sigma
rms_height_cm = sigma * 100.0
correlation_length_cm = corr_len * 100.0
soil_permittivity = complex(float(eps_r), float(eps_i))

k = 2 * np.pi / wavelength_m
ks = k * sigma
kl = k * corr_len

print(f"\nTest Case:")
print(f"  Frequency: {frequency_ghz} GHz")
print(f"  Wavelength: {wavelength_m*100:.2f} cm")
print(f"  Incidence: 40°")
print(f"  ℓ/σ ratio: {ratio:.1f}")
print(f"  εr: {eps_r:.1f} - j{eps_i:.1f}")
print(f"  σ: {sigma*100:.3f} cm (kσ = {ks:.3f})")
print(f"  ℓ: {corr_len*100:.3f} cm (kℓ = {kl:.3f})")

print(f"\nNMM3D Reference Values:")
print(f"  VV: {vv_ref:7.2f} dB")
print(f"  HH: {hh_ref:7.2f} dB")
print(f"  HV: {hv_ref:7.2f} dB")

# Compute AIEM values
print(f"\nAIEM Single Scattering (current implementation):")
for pol_name, pol_state in [("VV", PolarizationState.VV), 
                              ("HH", PolarizationState.HH), 
                              ("HV", PolarizationState.HV)]:
    result = mwRTMs.compute_soil_backscatter(
        model='aiem',
        radar_config=radar_config,
        frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm,
        correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity,
        correlation='exponential',
        polarization=pol_state,
    )
    result_db = 10.0 * np.log10(result) if result > 0 else -999
    print(f"  {pol_name}: {result_db:7.2f} dB (σ⁰ = {result:.6e})")

print("\n" + "=" * 80)
print("ANALYSIS:")
print("=" * 80)

# The issue: HV is essentially zero in AIEM single scattering
# This is EXPECTED because:
# 1. Single scattering from a rough surface with no tilted facets produces
#    negligible cross-polarization
# 2. Cross-pol arises primarily from:
#    a) Multiple scattering (what we're trying to add)
#    b) Tilted facets (geometric effect)
#    c) Volume scattering (not applicable for bare soil)

print("""
The catastrophic HV error (-305 dB bias) indicates that:

1. **AIEM single scattering produces near-zero HV** (as expected)
   - Single scattering from a smooth interface has no depolarization
   - The complementary term provides some HV but it's very small
   
2. **NMM3D includes multiple scattering effects** 
   - NMM3D is a full-wave numerical solution
   - It naturally includes all scattering orders
   - HV values in NMM3D are around -30 to -40 dB
   
3. **The multiple scattering module is not being called**
   - The test script doesn't enable multiple scattering by default
   - Need to use --add-multiple flag

SOLUTION: The multiple scattering implementation needs to be:
   a) Enabled in the test (--add-multiple flag)
   b) Validated to ensure it produces the correct HV contribution
   c) Checked for numerical accuracy (integration parameters)
""")

print("\nTesting with multiple scattering enabled:")
print("-" * 80)

try:
    hv_ms = mwRTMs.compute_soil_backscatter(
        model='aiem',
        radar_config=radar_config,
        frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm,
        correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity,
        correlation='exponential',
        polarization=PolarizationState.HV,
        include_multiple_scattering=True,
        ms_quadrature_points=65,  # Lower for speed
        ms_spectral_terms=6,
    )
    hv_ms_db = 10.0 * np.log10(hv_ms) if hv_ms > 0 else -999
    print(f"HV with MS: {hv_ms_db:7.2f} dB (σ⁰ = {hv_ms:.6e})")
    print(f"NMM3D HV:   {hv_ref:7.2f} dB")
    print(f"Difference: {hv_ms_db - hv_ref:7.2f} dB")
except Exception as e:
    print(f"Error computing MS: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 80)
