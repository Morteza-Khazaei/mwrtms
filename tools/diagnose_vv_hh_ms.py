"""Diagnostic script to compare VV and HH single vs multiple scattering."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry, PolarizationState
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Test parameters - from the error case
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
    enable_guardrails=False,  # Disable to see raw values
    ms_quadrature_points=129,
    ms_spectral_terms=8,
)

air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(er)

# Compute for both polarizations
print(f"\nTest case:")
print(f"  Frequency: {freq/1e9:.2f} GHz")
print(f"  Theta: {theta_deg}°")
print(f"  Sigma: {sigma*100:.2f} cm")
print(f"  L: {L*100:.2f} cm (L/sigma = {L/sigma:.1f})")
print(f"  ks = {wave.wavenumber * sigma:.3f}")
print(f"  er = {er}")

# Patch the model to extract single and MS contributions
import types

original_compute = model._compute_channel

def patched_compute(self, medium_above, medium_below, polarization, params):
    """Patched version that reports single and MS values."""
    # Call original but capture intermediate values
    # We'll do this by temporarily storing them
    self._last_single = None
    self._last_ms = None

    result = original_compute(medium_above, medium_below, polarization, params)
    return result

# Actually, let's just call the model twice - once with MS, once without
model_no_ms = AIEMModel(
    wave, geometry, surface,
    correlation_type="exponential",
    include_multiple_scattering=False,
    enable_guardrails=False,
)

print("\n" + "="*70)
print("VV POLARIZATION")
print("="*70)
vv_total = model.compute(air, soil, PolarizationState.VV)
vv_single = model_no_ms.compute(air, soil, PolarizationState.VV)
vv_ms = vv_total - vv_single

print(f"Single scattering:   {vv_single:.6e}  ({10*np.log10(vv_single):+.2f} dB)")
print(f"Multiple scattering: {vv_ms:.6e}  ({10*np.log10(vv_ms) if vv_ms > 0 else float('nan'):+.2f} dB)")
print(f"Total:               {vv_total:.6e}  ({10*np.log10(vv_total):+.2f} dB)")
print(f"MS/Single ratio:     {vv_ms/vv_single:.3f}  ({10*np.log10(vv_ms/vv_single):+.2f} dB)")

if vv_ms > vv_single:
    print(f"⚠️  WARNING: Multiple > Single by {vv_ms/vv_single:.2f}x")

print("\n" + "="*70)
print("HH POLARIZATION")
print("="*70)
hh_total = model.compute(air, soil, PolarizationState.HH)
hh_single = model_no_ms.compute(air, soil, PolarizationState.HH)
hh_ms = hh_total - hh_single

print(f"Single scattering:   {hh_single:.6e}  ({10*np.log10(hh_single):+.2f} dB)")
print(f"Multiple scattering: {hh_ms:.6e}  ({10*np.log10(hh_ms) if hh_ms > 0 else float('nan'):+.2f} dB)")
print(f"Total:               {hh_total:.6e}  ({10*np.log10(hh_total):+.2f} dB)")
print(f"MS/Single ratio:     {hh_ms/hh_single:.3f}  ({10*np.log10(hh_ms/hh_single):+.2f} dB)")

if hh_ms > hh_single:
    print(f"⚠️  WARNING: Multiple > Single by {hh_ms/hh_single:.2f}x")

print("\n" + "="*70)
print("COMPARISON")
print("="*70)
print(f"VV/HH single ratio:   {vv_single/hh_single:.3f}")
print(f"VV/HH multiple ratio: {vv_ms/hh_ms:.3f}")
print(f"VV MS is {vv_ms/hh_ms:.2f}x larger than HH MS")
