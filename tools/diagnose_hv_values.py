"""Diagnostic script to check HV values."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry, PolarizationState
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Test parameters
freq = 5.3e9  # Hz
theta_deg = 40.0
sigma = 0.015  # 1.5 cm RMS height
L = 0.06  # 6 cm correlation length (ratio = 4)
er = 12.0 + 1.8j

# Setup
wave = ElectromagneticWave(freq)
geometry = ScatteringGeometry(theta_i_deg=theta_deg)
surface = build_surface_from_statistics(sigma, L, correlation_type="exponential")

# Create model WITHOUT guardrails to get raw values
model = AIEMModel(
    wave, geometry, surface,
    correlation_type="exponential",
    include_multiple_scattering=True,
    enable_guardrails=False,
    ms_quadrature_points=129,
    ms_spectral_terms=8,
)

air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(er)

print(f"\nTest case:")
print(f"  Frequency: {freq/1e9:.2f} GHz")
print(f"  Theta: {theta_deg}°")
print(f"  Sigma: {sigma*100:.2f} cm")
print(f"  L: {L*100:.2f} cm (L/sigma = {L/sigma:.1f})")
print(f"  ks = {wave.wavenumber * sigma:.3f}")
print(f"  er = {er}")

print("\n" + "="*70)
print("HV POLARIZATION")
print("="*70)
hv_total = model.compute(air, soil, PolarizationState.HV)
print(f"HV (linear):     {hv_total:.6e}")
if hv_total > 0:
    print(f"HV (dB):         {10*np.log10(hv_total):+.2f} dB")
else:
    print(f"HV (dB):         NEGATIVE OR ZERO VALUE! {hv_total}")
