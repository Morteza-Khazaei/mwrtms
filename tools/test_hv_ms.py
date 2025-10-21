#!/usr/bin/env python3
"""Test HV multiple scattering specifically."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics

# Test parameters
frequency = 5.3e9  # Hz
theta_deg = 40.0
wavelength = 3e8 / frequency

# Moderate roughness
sigma_m = 0.01  # 1 cm RMS height
corr_length_m = 0.04  # 4 cm correlation length

# Create components
wave = ElectromagneticWave(frequency)
geometry = ScatteringGeometry(theta_i_deg=theta_deg)
surface = build_surface_from_statistics(
    rms_height_m=sigma_m,
    correlation_length_m=corr_length_m,
    correlation_type="exponential"
)

print(f"Testing HV multiple scattering:")
print(f"  Frequency: {frequency/1e9:.1f} GHz")
print(f"  Wavelength: {wavelength*100:.2f} cm")
print(f"  Incidence angle: {theta_deg}°")
print(f"  RMS height: {sigma_m*100:.1f} cm")
print(f"  Correlation length: {corr_length_m*100:.1f} cm")
print(f"  ks = {2*np.pi/wavelength * sigma_m:.3f}")
print(f"  kl = {2*np.pi/wavelength * corr_length_m:.3f}")
print()

# Test without MS
model_no_ms = AIEMModel(
    wave,
    geometry,
    surface,
    include_multiple_scattering=False,
    enable_guardrails=True
)

# Test with MS
model_with_ms = AIEMModel(
    wave,
    geometry,
    surface,
    include_multiple_scattering=True,
    ms_quadrature_points=65,
    ms_spectral_terms=6,
    enable_guardrails=True
)

air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(15.0 + 2.0j)

print("Without multiple scattering:")
result_no_ms = model_no_ms.run(air, soil)
print(f"  Result: {result_no_ms}")
print(f"  Type: {type(result_no_ms)}")
print(f"  Dir: {[x for x in dir(result_no_ms) if not x.startswith('_')]}")
print()

print("With multiple scattering:")
result_with_ms = model_with_ms.run(air, soil)
print(f"  Result: {result_with_ms}")
print()

# Try to extract values
try:
    data_no_ms = result_no_ms.as_dict()
    data_with_ms = result_with_ms.as_dict()
    
    print("Values without MS:")
    for pol in ['vv', 'hh', 'hv', 'vh']:
        val = data_no_ms.get(pol, 0.0)
        db = 10*np.log10(val) if val > 0 else -np.inf
        print(f"  {pol.upper()}: {db:.2f} dB ({val:.6e} linear)")
    
    print("\nValues with MS:")
    for pol in ['vv', 'hh', 'hv', 'vh']:
        val = data_with_ms.get(pol, 0.0)
        db = 10*np.log10(val) if val > 0 else -np.inf
        print(f"  {pol.upper()}: {db:.2f} dB ({val:.6e} linear)")
        
    print("\nMS contribution:")
    for pol in ['vv', 'hh', 'hv', 'vh']:
        val_no_ms = data_no_ms.get(pol, 0.0)
        val_with_ms = data_with_ms.get(pol, 0.0)
        ms_contrib = val_with_ms - val_no_ms
        pct = (ms_contrib / val_no_ms * 100) if val_no_ms > 0 else 0
        print(f"  {pol.upper()}: {ms_contrib:.6e} linear ({pct:+.1f}%)")
        
except Exception as e:
    print(f"Error extracting values: {e}")
    import traceback
    traceback.print_exc()
