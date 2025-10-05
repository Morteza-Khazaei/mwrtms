"""Test AIEM against a case with finite NMM3D HV value."""

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
k = 2 * np.pi / wavelength_m
radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)

# Find first case with finite HV and kσ > 0.5 (rougher surface)
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
kl = k * corr_len

print("=" * 80)
print("Testing AIEM HV with Finite NMM3D Reference")
print("=" * 80)

print(f"\nTest Parameters:")
print(f"  Frequency: {frequency_ghz} GHz (λ = {wavelength_m*100:.2f} cm)")
print(f"  Incidence: 40°")
print(f"  ℓ/σ: {ratio:.1f}")
print(f"  εr: {eps_r:.1f} - j{eps_i:.1f}")
print(f"  σ: {sigma*100:.3f} cm (kσ = {ks:.3f})")
print(f"  ℓ: {corr_len*100:.3f} cm (kℓ = {kl:.3f})")

print(f"\nNMM3D Reference:")
print(f"  VV: {vv_ref:7.2f} dB")
print(f"  HH: {hh_ref:7.2f} dB")
print(f"  HV: {hv_ref:7.2f} dB")

print(f"\nAIEM Single Scattering:")
vv_aiem = mwRTMs.compute_soil_backscatter(
    model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
    rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
    soil_permittivity=soil_permittivity, correlation='exponential',
    polarization=PolarizationState.VV,
)
hh_aiem = mwRTMs.compute_soil_backscatter(
    model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
    rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
    soil_permittivity=soil_permittivity, correlation='exponential',
    polarization=PolarizationState.HH,
)
hv_aiem = mwRTMs.compute_soil_backscatter(
    model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
    rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
    soil_permittivity=soil_permittivity, correlation='exponential',
    polarization=PolarizationState.HV,
)

vv_db = 10*np.log10(vv_aiem) if vv_aiem > 0 else -999
hh_db = 10*np.log10(hh_aiem) if hh_aiem > 0 else -999
hv_db = 10*np.log10(hv_aiem) if hv_aiem > 0 else -999

print(f"  VV: {vv_db:7.2f} dB (error: {vv_db-vv_ref:+.2f} dB)")
print(f"  HH: {hh_db:7.2f} dB (error: {hh_db-hh_ref:+.2f} dB)")
print(f"  HV: {hv_db:7.2f} dB (error: {hv_db-hv_ref:+.2f} dB) σ⁰={hv_aiem:.3e}")

print(f"\nAIEM with Multiple Scattering (65 pts, nmax=6):")
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
print(f"  HV: {hv_ms_db:7.2f} dB (error: {hv_ms_db-hv_ref:+.2f} dB) σ⁰={hv_ms:.3e}")

print(f"\nAIEM with Multiple Scattering (129 pts, nmax=8):")
hv_ms2 = mwRTMs.compute_soil_backscatter(
    model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
    rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
    soil_permittivity=soil_permittivity, correlation='exponential',
    polarization=PolarizationState.HV,
    include_multiple_scattering=True,
    ms_quadrature_points=129,
    ms_spectral_terms=8,
)
hv_ms2_db = 10*np.log10(hv_ms2) if hv_ms2 > 0 else -999
print(f"  HV: {hv_ms2_db:7.2f} dB (error: {hv_ms2_db-hv_ref:+.2f} dB) σ⁰={hv_ms2:.3e}")

print("\n" + "=" * 80)
print("CONCLUSION:")
print("=" * 80)
print("""
If multiple scattering still produces near-zero HV, the issues could be:

1. **Integration domain too small** - umax = 5/kℓ may not capture enough spectrum
2. **Exponential damping too strong** - exp(-σ²(...)) terms suppress contribution
3. **Wrong coordinate system** - backscatter geometry (φ_s = π vs φ_s = 0)
4. **Missing prefactors** - normalization constants incorrect
5. **Propagator signs** - sign errors in B coefficients or propagators
6. **Spectral function** - W_n(u,v,n) may have wrong normalization

The fact that VV and HH work well suggests the basic framework is correct,
but the cross-pol specific terms (B coefficients, cross-reflection R) may
have implementation errors.
""")
