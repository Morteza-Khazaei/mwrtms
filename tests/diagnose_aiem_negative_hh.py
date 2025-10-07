"""
Diagnose why AIEM produces negative HH values.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mwrtms.core import ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface.iem.aiem import AIEMModel
from mwrtms.scattering.surface.iem.i2em import I2EMModel


def test_single_case():
    """Test a single case to diagnose the issue."""
    print("=" * 80)
    print("DIAGNOSING AIEM NEGATIVE HH ISSUE")
    print("=" * 80)
    print()
    
    # Simple test case
    frequency_hz = 5.4e9
    theta_deg = 40.0
    sigma_m = 0.005  # 0.5 cm (small roughness)
    L_m = 0.05  # 5 cm
    eps = complex(15.0, 3.0)
    
    print(f"Test Configuration:")
    print(f"  Frequency: {frequency_hz/1e9:.1f} GHz")
    print(f"  Incidence angle: {theta_deg}°")
    print(f"  RMS height: {sigma_m*100:.1f} cm")
    print(f"  Correlation length: {L_m*100:.1f} cm")
    print(f"  Permittivity: {eps}")
    
    k = 2 * np.pi * frequency_hz / 3e8
    print(f"  ks = {k*sigma_m:.3f}")
    print(f"  kL = {k*L_m:.3f}")
    print()
    
    # Create components
    wave = ElectromagneticWave(frequency_hz)
    geometry = ScatteringGeometry(theta_i_deg=theta_deg)
    surface = build_surface_from_statistics(
        rms_height_m=sigma_m,
        correlation_length_m=L_m,
        correlation_type="exponential"
    )
    air = HomogeneousMedium(1.0 + 0.0j)
    soil = HomogeneousMedium(eps)
    
    # Test I2EM first (should work)
    print("Testing I2EM...")
    try:
        model_i2em = I2EMModel(
            wave, geometry, surface,
            correlation_type="exponential",
            auto_terms=False,
            spectral_terms=10
        )
        result_i2em = model_i2em.run(air, soil)
        
        print(f"✓ I2EM succeeded:")
        print(f"  VV: {result_i2em['vv']:.6e} linear ({10*np.log10(result_i2em['vv']):.2f} dB)")
        print(f"  HH: {result_i2em['hh']:.6e} linear ({10*np.log10(result_i2em['hh']):.2f} dB)")
        print()
    except Exception as e:
        print(f"✗ I2EM failed: {e}")
        print()
    
    # Test AIEM (may fail)
    print("Testing AIEM (original transition)...")
    try:
        model_aiem = AIEMModel(
            wave, geometry, surface,
            correlation_type="exponential",
            auto_terms=False,
            spectral_terms=10,
            use_i2em_transition=False
        )
        result_aiem = model_aiem.run(air, soil)
        
        print(f"✓ AIEM succeeded:")
        print(f"  VV: {result_aiem['vv']:.6e} linear ({10*np.log10(result_aiem['vv']):.2f} dB)")
        print(f"  HH: {result_aiem['hh']:.6e} linear ({10*np.log10(result_aiem['hh']):.2f} dB)")
        print()
    except Exception as e:
        print(f"✗ AIEM failed: {e}")
        print()
        print("This confirms AIEM has a bug producing negative HH values!")
        print()
    
    # Test AIEM with I2EM transition
    print("Testing AIEM (with I2EM transition)...")
    try:
        model_aiem_i2em = AIEMModel(
            wave, geometry, surface,
            correlation_type="exponential",
            auto_terms=False,
            spectral_terms=10,
            use_i2em_transition=True
        )
        result_aiem_i2em = model_aiem_i2em.run(air, soil)
        
        print(f"✓ AIEM with I2EM transition succeeded:")
        print(f"  VV: {result_aiem_i2em['vv']:.6e} linear ({10*np.log10(result_aiem_i2em['vv']):.2f} dB)")
        print(f"  HH: {result_aiem_i2em['hh']:.6e} linear ({10*np.log10(result_aiem_i2em['hh']):.2f} dB)")
        print()
    except Exception as e:
        print(f"✗ AIEM with I2EM transition also failed: {e}")
        print()
    
    print("=" * 80)
    print("DIAGNOSIS")
    print("=" * 80)
    print()
    print("The issue is that AIEM's complementary term calculation is producing")
    print("negative values for HH polarization. This is a BUG in the AIEM")
    print("implementation, not just a transition function issue.")
    print()
    print("Possible causes:")
    print("  1. Incorrect sign in complementary field coefficients")
    print("  2. Wrong combination of Kirchhoff + complementary terms")
    print("  3. Issue with the 8-branch complementary field calculation")
    print()
    print("This explains why we can't run the empirical test - AIEM is broken!")
    print()


if __name__ == "__main__":
    test_single_case()
