"""Diagnose multiple scattering balance issues."""

import sys
import math
import numpy as np

# Add source to path
sys.path.insert(0, '/home/morteza/usask/mwrtms/src')

from mwrtms.scattering.surface.iem.aiem import compute_aiem_backscatter

# Test parameters from NMM3D
frequency = 5.3e9  # Hz
wavelength = 3e8 / frequency  # meters
k = 2 * math.pi / wavelength

# Surface parameters for different ratios
test_cases = [
    {"ratio": 4, "sigma": 0.005, "corr_length": 0.020},   # ℓ/σ = 4
    {"ratio": 7, "sigma": 0.005, "corr_length": 0.035},   # ℓ/σ = 7
    {"ratio": 10, "sigma": 0.005, "corr_length": 0.050},  # ℓ/σ = 10
    {"ratio": 15, "sigma": 0.005, "corr_length": 0.075},  # ℓ/σ = 15
]

# Permittivity
eps_r = 15.0 + 3.5j

# Test angles
angles_deg = [20, 30, 40, 50, 60, 70]

print("=" * 100)
print("MULTIPLE SCATTERING BALANCE DIAGNOSTIC")
print("=" * 100)
print()

for case in test_cases:
    ratio = case["ratio"]
    sigma = case["sigma"]
    corr_length = case["corr_length"]
    ks = k * sigma
    
    print(f"\n{'='*100}")
    print(f"ℓ/σ = {ratio}, σ = {sigma*1000:.1f} mm, l = {corr_length*1000:.1f} mm, kσ = {ks:.3f}")
    print(f"{'='*100}")
    print(f"{'Angle':>6} {'VV_single':>12} {'VV_multi':>12} {'VV_ratio':>10} {'HH_single':>12} {'HH_multi':>12} {'HH_ratio':>10} {'HV_single':>12} {'HV_multi':>12} {'HV_ratio':>10}")
    print("-" * 100)
    
    violations = []
    
    for theta_deg in angles_deg:
        theta = math.radians(theta_deg)
        
        try:
            # Compute with multiple scattering
            result_with_ms = compute_aiem_backscatter(
                eps_r=eps_r,
                wavelength=wavelength,
                sigma_m=sigma,
                corr_length_m=corr_length,
                theta_i=theta,
                theta_s=theta,
                phi_i=0.0,
                phi_s=math.pi,
                surface_label="exponential",
                add_multiple_scattering=True,
            )
            
            # Compute without multiple scattering
            result_no_ms = compute_aiem_backscatter(
                eps_r=eps_r,
                wavelength=wavelength,
                sigma_m=sigma,
                corr_length_m=corr_length,
                theta_i=theta,
                theta_s=theta,
                phi_i=0.0,
                phi_s=math.pi,
                surface_label="exponential",
                add_multiple_scattering=False,
            )
            
            # Extract values
            vv_total = result_with_ms["vv"]
            vv_single = result_no_ms["vv"]
            vv_multi = vv_total - vv_single
            
            hh_total = result_with_ms["hh"]
            hh_single = result_no_ms["hh"]
            hh_multi = hh_total - hh_single
            
            hv_total = result_with_ms["hv"]
            hv_single = result_no_ms["hv"]
            hv_multi = hv_total - hv_single
            
            # Compute ratios
            vv_ratio = vv_multi / vv_single if vv_single > 0 else 0
            hh_ratio = hh_multi / hh_single if hh_single > 0 else 0
            hv_ratio = hv_multi / hv_single if abs(hv_single) > 1e-10 else float('inf')
            
            # Check for violations
            if vv_multi > vv_single:
                violations.append(f"  VV: θ={theta_deg}° MS/SS = {vv_ratio:.3f} > 1.0")
            if hh_multi > hh_single:
                violations.append(f"  HH: θ={theta_deg}° MS/SS = {hh_ratio:.3f} > 1.0")
            
            print(f"{theta_deg:>6} {vv_single:>12.6e} {vv_multi:>12.6e} {vv_ratio:>10.3f} "
                  f"{hh_single:>12.6e} {hh_multi:>12.6e} {hh_ratio:>10.3f} "
                  f"{hv_single:>12.6e} {hv_multi:>12.6e} {hv_ratio:>10.3f}")
            
        except Exception as e:
            print(f"{theta_deg:>6} ERROR: {str(e)[:80]}")
    
    if violations:
        print(f"\n⚠️  VIOLATIONS DETECTED for ℓ/σ = {ratio}:")
        for v in violations:
            print(v)
    else:
        print(f"\n✅ No violations for ℓ/σ = {ratio}")

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
print()
print("Rule C.1 (Yang et al. 2017): For co-pol, multiple scattering should NOT exceed single scattering")
print("  - For kσ < 0.3: Single >> Multiple")
print("  - For kσ ~ 0.5: Single > Multiple (but comparable)")
print("  - For kσ > 1.0: Single ~ Multiple")
print()
print("If violations occur, possible causes:")
print("  1. Integration domain too large (including evanescent modes)")
print("  2. Spectral summation not converged")
print("  3. Field coefficients incorrect")
print("  4. Phase factors wrong sign")
print("  5. Normalization error")
