"""Quick test to check multiple scattering magnitude."""
import numpy as np
from pathlib import Path
from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState

# Load one test case
lut_path = Path("data/NMM3D_LUT_NRCS_40degree.dat")
table = np.loadtxt(lut_path)
mask = np.isclose(table[:, 0], 40.0)
rows = table[mask]

# Find case with finite HV and kσ > 0.5
frequency_ghz = 5.405
wavelength_m = 0.3 / frequency_ghz
k = 2 * np.pi / wavelength_m
radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)

for row in rows:
    theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
    ks = k * rms_norm * wavelength_m
    if np.isfinite(hv_ref) and ks > 0.5:
        break

sigma = rms_norm * wavelength_m
corr_len = ratio * sigma
rms_height_cm = sigma * 100.0
correlation_length_cm = corr_len * 100.0
soil_permittivity = complex(float(eps_r), float(eps_i))

print(f"Test case: kσ={ks:.3f}, kℓ={k*corr_len:.3f}, ℓ/σ={ratio:.1f}")
print(f"NMM3D HV: {hv_ref:.2f} dB")

# Test with MS
hv_ms = mwRTMs.compute_soil_backscatter(
    model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
    rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
    soil_permittivity=soil_permittivity, correlation='exponential',
    polarization=PolarizationState.HV,
    include_multiple_scattering=True,
    ms_quadrature_points=65,
    ms_spectral_terms=6,
)

hv_ms_db = 10*np.log10(hv_ms) if hv_ms > 0 else -999
error_db = hv_ms_db - hv_ref
error_linear_factor = 10**(error_db/10)

print(f"AIEM HV with MS: {hv_ms_db:.2f} dB (σ⁰={hv_ms:.3e})")
print(f"Error: {error_db:.2f} dB")
print(f"Linear factor off by: {error_linear_factor:.3e} (need to multiply by {1/error_linear_factor:.3e})")
print(f"In dB, need to add: {-error_db:.2f} dB")

# Check what factor we need
needed_factor = 10**(-error_db/20)  # sqrt because sigma0 ~ |field|^2
print(f"\nField amplitude needs factor of: {needed_factor:.3e}")
print(f"Power needs factor of: {needed_factor**2:.3e}")
