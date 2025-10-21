"""Test cross-pol guardrails are working."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry, PolarizationState
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Test parameters
freq = 5.3e9  # Hz
theta_deg = 40.0
sigma = 0.015  # 1.5 cm RMS height
L = 0.06  # 6 cm correlation length
er = 12.0 + 1.8j

# Setup
wave = ElectromagneticWave(freq)
geometry = ScatteringGeometry(theta_i_deg=theta_deg)
surface = build_surface_from_statistics(sigma, L, correlation_type="exponential")

# Create model WITH guardrails
model = AIEMModel(
    wave, geometry, surface,
    correlation_type="exponential",
    include_multiple_scattering=True,
    enable_guardrails=True,  # Enable guardrails
    ms_quadrature_points=129,
    ms_spectral_terms=8,
)

air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(er)

print(f"\nTesting Cross-Pol Guardrails")
print(f"="*70)

print(f"\nTest case:")
print(f"  Frequency: {freq/1e9:.2f} GHz")
print(f"  Theta: {theta_deg}°")
print(f"  Sigma: {sigma*100:.2f} cm")
print(f"  L: {L*100:.2f} cm")
print(f"  ks = {wave.wavenumber * sigma:.3f}")

print(f"\nComputing HV with guardrails enabled...")
try:
    hv_total = model.compute(air, soil, PolarizationState.HV)
    hv_db = 10*np.log10(hv_total) if hv_total > 0 else float('-inf')
    print(f"✓ HV computation successful")
    print(f"  HV (linear): {hv_total:.6e}")
    print(f"  HV (dB):     {hv_db:.2f} dB")

    if hv_total < 1e-10:
        print(f"  ⚠️  Very small value - check if physical")
    elif hv_db > 0:
        print(f"  ⚠️  Above 0 dB - unusually large for natural surface")
    else:
        print(f"  ✓ Value in reasonable range")

except Exception as e:
    print(f"✗ Guardrail violation detected:")
    print(f"  {type(e).__name__}: {e}")

print(f"\n" + "="*70)
print("Guardrail test complete!")
