#!/usr/bin/env python3
"""Debug HV multiple scattering values."""

import numpy as np
import math

# Import the multiple scattering function directly
from src.mwrtms.scattering.iem.multiple_scattering import compute_multiple_scattering

# Test parameters from NMM3D (typical case)
frequency_ghz = 5.405
wavelength = 3e8 / (frequency_ghz * 1e9)  # meters
k = 2 * np.pi / wavelength

# Surface parameters
rms_height_cm = 0.5
correlation_length_cm = 2.0  # ratio = 4
sigma = rms_height_cm / 100.0  # meters
L = correlation_length_cm / 100.0  # meters

ks = k * sigma
kl = k * L

# Geometry
theta_i = np.radians(40.0)
theta_s = theta_i  # backscatter
phi_i = 0.0
phi_s = np.pi  # backscatter

# Soil permittivity
eps_r = complex(3.0, 1.0)

print("=" * 60)
print("HV Multiple Scattering Diagnostic")
print("=" * 60)
print(f"Frequency: {frequency_ghz} GHz")
print(f"Wavelength: {wavelength*100:.4f} cm")
print(f"k: {k:.2f} rad/m")
print(f"RMS height: {rms_height_cm} cm (σ = {sigma*100:.4f} cm)")
print(f"Correlation length: {correlation_length_cm} cm (L = {L*100:.4f} cm)")
print(f"ks: {ks:.4f}")
print(f"kl: {kl:.4f}")
print(f"Ratio l/σ: {L/sigma:.1f}")
print(f"Incidence angle: 40°")
print(f"Permittivity: {eps_r}")
print()

# Compute multiple scattering for all polarizations
print("Computing multiple scattering...")
print()

for pol in ['vv', 'hh', 'hv']:
    try:
        ms_result = compute_multiple_scattering(
            theta_i=theta_i,
            theta_s=theta_s,
            phi_i=phi_i,
            phi_s=phi_s,
            er=eps_r,
            ks=ks,
            kl=kl,
            k=k,
            sigma=sigma,
            corr_length=L,
            surface_label='exponential',
            polarisations=[pol],
            n_points=129,
            nmax=8
        )
        
        sigma_ms = ms_result[pol]
        
        if sigma_ms > 0:
            sigma_ms_db = 10 * np.log10(sigma_ms)
        else:
            sigma_ms_db = float('-inf')
        
        print(f"{pol.upper():3s}: σ_ms = {sigma_ms:.6e} (linear)")
        print(f"     σ_ms = {sigma_ms_db:+.2f} dB")
        print()
        
    except Exception as e:
        print(f"{pol.upper():3s}: ERROR - {e}")
        import traceback
        traceback.print_exc()
        print()

print("=" * 60)
print("Analysis:")
print("=" * 60)
print("If HV shows -inf dB or very small values (< 1e-20),")
print("then multiple scattering is essentially zero for cross-pol.")
print()
print("Expected behavior:")
print("- VV/HH: Should have finite MS contribution (typically -30 to -50 dB)")
print("- HV: Should also have finite MS contribution (typically -40 to -60 dB)")
print()
print("If HV is near zero, check:")
print("1. Are the propagators being computed correctly for cross-pol?")
print("2. Are the B coefficients correct?")
print("3. Is the integration producing non-zero results?")
