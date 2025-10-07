"""
Debug the transition function to see what's happening.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mwrtms.scattering.surface.iem.transition import compute_transition_function
from mwrtms.scattering.surface.iem.transition_i2em import compute_i2em_transition_function
from mwrtms.scattering.surface.iem.spectrum_aiem import compute_aiem_spectrum
from mwrtms.scattering.surface.iem.geometry_utils import compute_spatial_frequency


def test_transition_functions():
    """Compare the two transition functions."""
    print("=" * 80)
    print("Transition Function Comparison")
    print("=" * 80)
    print()
    
    # Test parameters
    eps_r = complex(15.0, 3.0)
    theta_i = np.radians(40.0)
    k = 2 * np.pi * 5.4e9 / 3e8  # wavenumber
    sigma = 0.01  # RMS height
    L = 0.10  # correlation length
    ks = k * sigma
    cs = np.cos(theta_i)
    
    print(f"Parameters:")
    print(f"  eps_r = {eps_r}")
    print(f"  theta_i = {np.degrees(theta_i):.1f}°")
    print(f"  ks = {ks:.3f}")
    print(f"  kL = {k*L:.3f}")
    print()
    
    # Compute spectrum
    n_terms = 15
    kl = k * L
    K = compute_spatial_frequency(kl, theta_i, theta_i, np.pi, 0.0)
    
    spectra = np.zeros(n_terms)
    for n in range(1, n_terms + 1):
        spectra[n-1] = compute_aiem_spectrum(kl, K, n, "exponential", 1.5)
    
    print(f"Computed {n_terms} spectral terms")
    print(f"  First 5 spectrum values: {spectra[:5]}")
    print()
    
    # Test AIEM transition
    print("AIEM Transition Function:")
    print("-" * 80)
    Tfv_aiem, Tfh_aiem = compute_transition_function(
        eps_r, theta_i, ks, cs, spectra, n_terms
    )
    print(f"  Tfv = {Tfv_aiem:.6f}")
    print(f"  Tfh = {Tfh_aiem:.6f}")
    print()
    
    # Test I2EM transition
    print("I2EM Transition Function:")
    print("-" * 80)
    Tfv_i2em, Tfh_i2em = compute_i2em_transition_function(
        eps_r, theta_i, ks, cs, spectra, n_terms
    )
    print(f"  Tfv = {Tfv_i2em:.6f}")
    print(f"  Tfh = {Tfh_i2em:.6f}")
    print()
    
    # Compare
    print("=" * 80)
    print("Comparison:")
    print("=" * 80)
    print(f"  ΔTfv = {Tfv_aiem - Tfv_i2em:+.6f}")
    print(f"  ΔTfh = {Tfh_aiem - Tfh_i2em:+.6f}")
    print()
    
    # Test Fresnel coefficients
    from mwrtms.scattering.surface.iem.fresnel_utils import (
        compute_fresnel_incident,
        compute_fresnel_nadir
    )
    
    Rvi, Rhi, Rvhi = compute_fresnel_incident(eps_r, theta_i)
    rv0, rh0 = compute_fresnel_nadir(eps_r)
    
    print("Fresnel Coefficients:")
    print("-" * 80)
    print(f"  Rvi = {Rvi:.6f}")
    print(f"  Rhi = {Rhi:.6f}")
    print(f"  rv0 = {rv0:.6f}")
    print(f"  rh0 = {rh0:.6f}")
    print()
    
    # Compute transition-adjusted coefficients
    print("Transition-Adjusted Reflection Coefficients:")
    print("-" * 80)
    
    # AIEM style
    Rvtran_aiem = Rvi + (rv0 - Rvi) * Tfv_aiem
    Rhtran_aiem = Rhi + (rh0 - Rhi) * Tfh_aiem
    print(f"AIEM:")
    print(f"  Rvtran = {Rvtran_aiem:.6f}")
    print(f"  Rhtran = {Rhtran_aiem:.6f}")
    print()
    
    # I2EM style (uses negative rh0)
    rh0_i2em = -rv0
    Rvtran_i2em = Rvi + (rv0 - Rvi) * Tfv_i2em
    Rhtran_i2em = Rhi + (rh0_i2em - Rhi) * Tfh_i2em
    print(f"I2EM Transition (rh0 = {rh0_i2em:.6f}):")
    print(f"  Rvtran = {Rvtran_i2em:.6f}")
    print(f"  Rhtran = {Rhtran_i2em:.6f}")
    print()
    
    print("=" * 80)
    print("ANALYSIS:")
    print("=" * 80)
    
    if Tfv_i2em == Tfh_i2em:
        print("✓ I2EM uses same transition for both polarizations")
    else:
        print("✗ I2EM transitions differ (unexpected!)")
    
    if abs(Tfv_aiem - Tfh_aiem) > 0.01:
        print("✓ AIEM uses different transitions for V and H")
    else:
        print("✗ AIEM transitions are similar")
    
    print()


if __name__ == "__main__":
    test_transition_functions()
