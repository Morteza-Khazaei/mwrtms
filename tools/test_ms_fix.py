#!/usr/bin/env python3
"""Quick test to verify multiple scattering fix."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics

# Test parameters from the failing case
frequency = 5.3e9  # Hz
theta_deg = 40.0
wavelength = 3e8 / frequency

# High roughness case that was failing
sigma_m = 0.02  # 2 cm RMS height
corr_length_m = 0.08  # 8 cm correlation length (ratio=4.0)

# Create components
wave = ElectromagneticWave(frequency)
geometry = ScatteringGeometry(theta_i_deg=theta_deg)
surface = build_surface_from_statistics(
    rms_height_m=sigma_m,
    correlation_length_m=corr_length_m,
    correlation_type="exponential"
)

# Create model with multiple scattering enabled
model = AIEMModel(
    wave,
    geometry,
    surface,
    include_multiple_scattering=True,
    ms_quadrature_points=65,  # Reduced for speed
    ms_spectral_terms=6,
    enable_guardrails=True
)

# Media
air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(15.0 + 2.0j)

print(f"Testing AIEM with multiple scattering:")
print(f"  Frequency: {frequency/1e9:.1f} GHz")
print(f"  Wavelength: {wavelength*100:.2f} cm")
print(f"  Incidence angle: {theta_deg}°")
print(f"  RMS height: {sigma_m*100:.1f} cm")
print(f"  Correlation length: {corr_length_m*100:.1f} cm")
print(f"  ks = {2*np.pi/wavelength * sigma_m:.3f}")
print(f"  kl = {2*np.pi/wavelength * corr_length_m:.3f}")
print()

try:
    result = model.run(air, soil)
    print("✓ Computation successful!")
    print(f"Result type: {type(result)}")
    print(f"Result: {result}")
        
except Exception as e:
    print(f"✗ Computation failed: {e}")
    import traceback
    traceback.print_exc()
