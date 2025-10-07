"""
Diagnose why HH gets worse with I2EM transition.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mwrtms.scattering.surface.iem.fresnel_utils import (
    compute_fresnel_incident,
    compute_fresnel_nadir
)


def analyze_hh_problem():
    """Analyze what happens to HH reflection coefficients."""
    print("=" * 80)
    print("DIAGNOSING HH PROBLEM WITH I2EM TRANSITION")
    print("=" * 80)
    print()
    
    # Test case
    eps_r = complex(15.0, 3.0)
    theta_i = np.radians(40.0)
    
    # Compute Fresnel coefficients
    Rvi, Rhi, Rvhi = compute_fresnel_incident(eps_r, theta_i)
    rv0, rh0 = compute_fresnel_nadir(eps_r)
    
    print(f"Fresnel Coefficients:")
    print(f"  Rvi = {Rvi:.6f}")
    print(f"  Rhi = {Rhi:.6f}")
    print(f"  rv0 = {rv0:.6f}")
    print(f"  rh0 = {rh0:.6f} (AIEM convention)")
    print()
    
    # Simulate transition factors
    Tf_aiem_v = 0.81
    Tf_aiem_h = 0.93
    Tf_i2em = 0.35
    
    print("=" * 80)
    print("AIEM ORIGINAL TRANSITION")
    print("=" * 80)
    print()
    print(f"Transition factors: Tfv = {Tf_aiem_v}, Tfh = {Tf_aiem_h}")
    print(f"Using rh0 = {rh0:.6f} (positive)")
    print()
    
    Rvtran_aiem = Rvi + (rv0 - Rvi) * Tf_aiem_v
    Rhtran_aiem = Rhi + (rh0 - Rhi) * Tf_aiem_h
    
    print(f"VV: Rvtran = {Rvtran_aiem:.6f}")
    print(f"HH: Rhtran = {Rhtran_aiem:.6f}")
    print()
    
    print("=" * 80)
    print("I2EM TRANSITION (CURRENT IMPLEMENTATION)")
    print("=" * 80)
    print()
    print(f"Transition factor: Tf = {Tf_i2em} (same for both)")
    print(f"Using rh0_trans = -rv0 = {-rv0:.6f} (NEGATIVE!)")
    print()
    
    rh0_i2em = -rv0
    Rvtran_i2em = Rvi + (rv0 - Rvi) * Tf_i2em
    Rhtran_i2em = Rhi + (rh0_i2em - Rhi) * Tf_i2em
    
    print(f"VV: Rvtran = {Rvtran_i2em:.6f}")
    print(f"HH: Rhtran = {Rhtran_i2em:.6f}")
    print()
    
    print("=" * 80)
    print("THE PROBLEM")
    print("=" * 80)
    print()
    
    print(f"HH Reflection Coefficient Change:")
    print(f"  AIEM:  Rhi = {Rhi:.6f} → Rhtran = {Rhtran_aiem:.6f}")
    print(f"  I2EM:  Rhi = {Rhi:.6f} → Rhtran = {Rhtran_i2em:.6f}")
    print()
    
    print(f"The issue:")
    print(f"  • Rhi is NEGATIVE ({Rhi:.3f})")
    print(f"  • With AIEM: transitions toward POSITIVE rh0 ({rh0:.3f})")
    print(f"    → Rhtran becomes positive ({Rhtran_aiem:.3f})")
    print()
    print(f"  • With I2EM: transitions toward NEGATIVE rh0 ({rh0_i2em:.3f})")
    print(f"    → Rhtran stays VERY negative ({Rhtran_i2em:.3f})")
    print()
    
    print(f"Impact on backscatter:")
    print(f"  fhh ∝ Rhtran")
    print(f"  σ_hh ∝ |fhh|²")
    print()
    print(f"  AIEM: |Rhtran|² = {abs(Rhtran_aiem)**2:.6f}")
    print(f"  I2EM: |Rhtran|² = {abs(Rhtran_i2em)**2:.6f}")
    print()
    print(f"  Ratio: {abs(Rhtran_i2em)**2 / abs(Rhtran_aiem)**2:.3f}")
    print(f"  In dB: {10*np.log10(abs(Rhtran_i2em)**2 / abs(Rhtran_aiem)**2):.2f} dB")
    print()
    
    print("=" * 80)
    print("ROOT CAUSE")
    print("=" * 80)
    print()
    print("The I2EM transition uses rh0 = -rv0, which is correct for I2EM.")
    print("But AIEM's complementary terms expect rh0 = +rv0!")
    print()
    print("When we use I2EM's negative rh0 in AIEM:")
    print("  1. Rhtran becomes much more negative")
    print("  2. This affects fhh (Kirchhoff coefficient)")
    print("  3. The complementary terms still use Rhi (not Rhtran)")
    print("  4. There's a mismatch between Kirchhoff and complementary terms")
    print()
    print("This causes HH to get WORSE instead of better!")
    print()
    
    print("=" * 80)
    print("THE FIX")
    print("=" * 80)
    print()
    print("Option 1: Don't use negative rh0 for HH in AIEM")
    print("  → Use I2EM transition factor but keep rh0 = +rv0")
    print()
    print("Option 2: Only use I2EM transition for VV, keep AIEM for HH")
    print("  → Hybrid approach")
    print()
    print("Option 3: Fix the complementary terms to work with negative rh0")
    print("  → More complex, requires understanding complementary term physics")
    print()
    print("RECOMMENDATION: Option 1 - Use I2EM Tf but keep rh0 positive")
    print()


if __name__ == "__main__":
    analyze_hh_problem()
