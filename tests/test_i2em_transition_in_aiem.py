"""
Test AIEM with I2EM transition method vs original AIEM transition.
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
    """Test a single representative case."""
    print("=" * 80)
    print("Testing AIEM with I2EM Transition Method")
    print("=" * 80)
    print()
    
    # Test parameters (typical soil conditions)
    frequency_hz = 5.4e9  # C-band
    theta_deg = 40.0
    sigma_m = 0.01  # 1 cm RMS height
    L_m = 0.10  # 10 cm correlation length
    eps = complex(15.0, 3.0)  # Typical soil permittivity
    
    print(f"Test Configuration:")
    print(f"  Frequency: {frequency_hz/1e9:.1f} GHz")
    print(f"  Incidence angle: {theta_deg}°")
    print(f"  RMS height: {sigma_m*100:.1f} cm")
    print(f"  Correlation length: {L_m*100:.1f} cm")
    print(f"  Permittivity: {eps.real:.1f} + {eps.imag:.1f}j")
    print()
    
    # Create wave and geometry
    wave = ElectromagneticWave(frequency_hz)
    geometry = ScatteringGeometry(theta_i_deg=theta_deg)
    
    # Create surface
    surface = build_surface_from_statistics(
        rms_height_m=sigma_m,
        correlation_length_m=L_m,
        correlation_type="exponential"
    )
    
    # Create media
    air = HomogeneousMedium(1.0 + 0.0j)
    soil = HomogeneousMedium(eps)
    
    print("Computing backscatter coefficients...")
    print()
    
    # Test 1: AIEM with original transition
    print("1. AIEM with Original Transition:")
    model_orig = AIEMModel(
        wave, geometry, surface,
        correlation_type="exponential",
        auto_terms=False,
        spectral_terms=15,
        use_i2em_transition=False
    )
    result_orig = model_orig.run(air, soil)
    result_orig_db = result_orig.to_db()
    
    print(f"   VV: {result_orig_db['vv']:7.2f} dB")
    print(f"   HH: {result_orig_db['hh']:7.2f} dB")
    print(f"   HV: {result_orig_db['hv']:7.2f} dB")
    print()
    
    # Test 2: AIEM with I2EM transition
    print("2. AIEM with I2EM Transition:")
    model_i2em_trans = AIEMModel(
        wave, geometry, surface,
        correlation_type="exponential",
        auto_terms=False,
        spectral_terms=15,
        use_i2em_transition=True
    )
    result_i2em_trans = model_i2em_trans.run(air, soil)
    result_i2em_trans_db = result_i2em_trans.to_db()
    
    print(f"   VV: {result_i2em_trans_db['vv']:7.2f} dB")
    print(f"   HH: {result_i2em_trans_db['hh']:7.2f} dB")
    print(f"   HV: {result_i2em_trans_db['hv']:7.2f} dB")
    print()
    
    # Test 3: I2EM model for reference
    print("3. I2EM Model (Reference):")
    model_i2em = I2EMModel(
        wave, geometry, surface,
        correlation_type="exponential",
        auto_terms=False,
        spectral_terms=15
    )
    result_i2em = model_i2em.run(air, soil)
    result_i2em_db = result_i2em.to_db()
    
    print(f"   VV: {result_i2em_db['vv']:7.2f} dB")
    print(f"   HH: {result_i2em_db['hh']:7.2f} dB")
    print(f"   HV: {result_i2em_db['hv']:7.2f} dB")
    print()
    
    # Compute differences
    print("=" * 80)
    print("Differences (AIEM Original - AIEM I2EM Transition):")
    print("=" * 80)
    diff_vv = result_orig_db['vv'] - result_i2em_trans_db['vv']
    diff_hh = result_orig_db['hh'] - result_i2em_trans_db['hh']
    diff_hv = result_orig_db['hv'] - result_i2em_trans_db['hv']
    
    print(f"   VV: {diff_vv:+7.2f} dB")
    print(f"   HH: {diff_hh:+7.2f} dB")
    print(f"   HV: {diff_hv:+7.2f} dB")
    print()
    
    print("Differences (AIEM I2EM Transition - I2EM Model):")
    print("-" * 80)
    diff_vv_i2em = result_i2em_trans_db['vv'] - result_i2em_db['vv']
    diff_hh_i2em = result_i2em_trans_db['hh'] - result_i2em_db['hh']
    diff_hv_i2em = result_i2em_trans_db['hv'] - result_i2em_db['hv']
    
    print(f"   VV: {diff_vv_i2em:+7.2f} dB")
    print(f"   HH: {diff_hh_i2em:+7.2f} dB")
    print(f"   HV: {diff_hv_i2em:+7.2f} dB")
    print()
    
    print("=" * 80)
    print("ANALYSIS:")
    print("=" * 80)
    
    if abs(diff_vv) > 0.5 or abs(diff_hh) > 0.5:
        print("✓ I2EM transition method CHANGES AIEM co-pol results significantly!")
        print(f"  - VV changed by {diff_vv:+.2f} dB")
        print(f"  - HH changed by {diff_hh:+.2f} dB")
    else:
        print("✗ I2EM transition method has minimal effect on AIEM co-pol results.")
    
    print()
    
    if abs(diff_vv_i2em) < 1.0 and abs(diff_hh_i2em) < 1.0:
        print("✓ AIEM with I2EM transition is CLOSE to I2EM model!")
        print(f"  - VV difference: {diff_vv_i2em:+.2f} dB")
        print(f"  - HH difference: {diff_hh_i2em:+.2f} dB")
    else:
        print("✗ AIEM with I2EM transition still differs from I2EM model.")
        print(f"  - VV difference: {diff_vv_i2em:+.2f} dB")
        print(f"  - HH difference: {diff_hh_i2em:+.2f} dB")
        print("  (This is expected due to complementary term differences)")
    
    print()


if __name__ == "__main__":
    test_single_case()
