"""
Diagnostic tool to check sign convention consistency between VKA and complementary terms.

This compares the wave vector definitions used in:
1. VKA (kirchhoff.py) - the reference implementation
2. Complementary terms (complementary.py)
"""

import numpy as np
import sys
sys.path.insert(0, 'src')

from mwrtms.scattering.surface.iem.kirchhoff import VKA

def analyze_vka_convention():
    """Analyze the VKA wave vector convention for backscatter."""

    # Backscatter geometry at 40 degrees
    theta_i = np.deg2rad(40.0)
    theta_s = np.deg2rad(40.0)
    phi_i = 0.0
    phi_s = np.pi  # Backscatter condition

    # Dummy Fresnel coefficients
    Rv = 0.5 + 0.1j
    Rh = 0.4 + 0.1j

    # Create VKA instance
    vka = VKA(theta_i, theta_s, phi_i, phi_s, Rv, Rh)

    # Get all vectors
    vectors = vka.get_all_vectors()

    print("="*70)
    print("VKA WAVE VECTOR CONVENTION (REFERENCE)")
    print("="*70)
    print("\nBackscatter geometry: θ_i = θ_s = 40°, φ_s = 180°")
    print("\n--- Incident wave vector k_i ---")
    print(f"k_i = [{vectors['k_i'][0]:+.6f}, {vectors['k_i'][1]:+.6f}, {vectors['k_i'][2]:+.6f}]")
    print(f"  sin(θ_i) = {np.sin(theta_i):.6f}")
    print(f"  cos(θ_i) = {np.cos(theta_i):.6f}")
    print(f"  k_i[2] = -cos(θ_i) = {-np.cos(theta_i):.6f}  ← DOWNWARD propagation (negative z)")

    print("\n--- Scattered wave vector k_s ---")
    print(f"k_s = [{vectors['k_s'][0]:+.6f}, {vectors['k_s'][1]:+.6f}, {vectors['k_s'][2]:+.6f}]")
    print(f"  sin(θ_s) = {np.sin(theta_s):.6f}")
    print(f"  cos(θ_s) = {np.cos(theta_s):.6f}")
    print(f"  k_s[2] = +cos(θ_s) = {+np.cos(theta_s):.6f}  ← UPWARD propagation (positive z)")

    print("\n--- Surface slopes (stationary phase) ---")
    zx, zy = vectors['surface_slopes']
    print(f"zx = -(k_sx - k_ix)/(k_sz - k_iz) = {zx:.6f}")
    print(f"zy = -(k_sy - k_iy)/(k_sz - k_iz) = {zy:.6f}")

    kx_i = vectors['k_i'][0]
    ky_i = vectors['k_i'][1]
    kz_i = vectors['k_i'][2]
    kx_s = vectors['k_s'][0]
    ky_s = vectors['k_s'][1]
    kz_s = vectors['k_s'][2]

    print(f"\nVerification:")
    print(f"  k_sx - k_ix = {kx_s - kx_i:.6f}")
    print(f"  k_sy - k_iy = {ky_s - ky_i:.6f}")
    print(f"  k_sz - k_iz = {kz_s - kz_i:.6f}")
    print(f"  zx = -({kx_s - kx_i:.6f})/({kz_s - kz_i:.6f}) = {-(kx_s - kx_i)/(kz_s - kz_i):.6f}")
    print(f"  zy = -({ky_s - ky_i:.6f})/({kz_s - kz_i:.6f}) = {-(ky_s - ky_i)/(kz_s - kz_i):.6f}")

    print("\n" + "="*70)
    print("CRITICAL SIGN CONVENTION")
    print("="*70)
    print(f"k_i[2] = {kz_i:+.6f}  (NEGATIVE: downward)")
    print(f"k_s[2] = {kz_s:+.6f}  (POSITIVE: upward)")
    print(f"\nFor complementary terms with q (vertical wave component):")
    print(f"  Air-side upward:   +q  (matches k_s sign)")
    print(f"  Air-side downward: -q  (matches k_i sign)")

    return vectors


def check_complementary_conventions():
    """Check what convention is being used in complementary terms."""

    print("\n" + "="*70)
    print("COMPLEMENTARY TERMS CONVENTION (aiem.py)")
    print("="*70)

    print("\nFrom aiem.py line ~504-565, the complementary terms are called with:")
    print("\nAir-side incident upward (Fvaupi):")
    print("  u = -si  (negative sin(θ_i))")
    print("  q = qq1  (positive)")
    print("  qslp = qq1  (positive)")
    print("  direction = +1")
    print("  → Should match k_s (upward, positive z)")

    print("\nAir-side incident downward (Fvadni):")
    print("  u = -si  (negative sin(θ_i))")
    print("  q = -qq1  (NEGATIVE)")
    print("  qslp = -qq1  (NEGATIVE)")
    print("  direction = -1")
    print("  → Should match k_i (downward, negative z)")

    print("\nAir-side scattered upward (Fvaups):")
    print("  u = -sis*csfs  (for backscatter: -sin(θ_s)*cos(π))")
    print("  v = -sis*sfs   (for backscatter: -sin(θ_s)*sin(π))")
    print("  q = qq2  (positive)")
    print("  qslp = qq2  (positive)")
    print("  direction = +1")
    print("  → Should match k_s (upward, positive z)")

    print("\nAir-side scattered downward (Fvadns):")
    print("  u = -sis*csfs")
    print("  v = -sis*sfs")
    print("  q = -qq2  (NEGATIVE)")
    print("  qslp = -qq2  (NEGATIVE)")
    print("  direction = -1")
    print("  → Should match k_i (downward, negative z)")

    print("\n" + "="*70)
    print("SIGN ISSUE DIAGNOSIS")
    print("="*70)

    print("\nThe 'direction' parameter in _signed_qfix():")
    print("  direction = +1  →  returns +qfix  (upward)")
    print("  direction = -1  →  returns -qfix  (downward)")

    print("\nIn complementary.py, the qfix denominator is used as:")
    print("  av = (1 + R) / signed_qfix")
    print("  bv = (1 - R) / signed_qfix")

    print("\nFOR BACKSCATTER:")
    print("  φ_s = φ_i + π  →  cos(φ_s) = -cos(φ_i), sin(φ_s) = -sin(φ_i)")
    print("  For φ_i = 0:")
    print("    cos(φ_s) = -1")
    print("    sin(φ_s) = 0")
    print("    sis*csfs = sin(θ_s)*(-1) = -sin(θ_s)")
    print("    sis*sfs = sin(θ_s)*(0) = 0")


def main():
    print("\n" + "="*70)
    print("AIEM SIGN CONVENTION DIAGNOSTIC")
    print("="*70)

    # Analyze VKA convention
    vka_vectors = analyze_vka_convention()

    # Check complementary conventions
    check_complementary_conventions()

    print("\n" + "="*70)
    print("RECOMMENDATION")
    print("="*70)
    print("\n1. VKA defines:")
    print("   - k_i[2] = -cos(θ_i)  (downward)")
    print("   - k_s[2] = +cos(θ_s)  (upward)")

    print("\n2. Complementary terms should follow the SAME convention:")
    print("   - Upward propagation: positive q, direction=+1")
    print("   - Downward propagation: negative q, direction=-1")

    print("\n3. CHECK: Are the slope calculations consistent?")
    print("   In complementary.py _compute_slope_and_field_components():")
    print("   - zx = -ksxu / (css - qslp)")
    print("   - zy = -ksyv / (css - qslp)")
    print("\n   These should match VKA stationary phase slopes!")

    print("\n4. CRITICAL: The cterm (c1-c6) formulas must be derived with")
    print("   the SAME sign convention as VKA!")
    print("\n   Currently seeing +4.7 dB overestimation in VV")
    print("   and +8.7 dB overestimation in HH.")
    print("\n   This suggests a factor of ~2-3x error, which could be:")
    print("   - Wrong sign in a term that should cancel")
    print("   - Missing negative sign in reflection coefficient combination")
    print("   - Inconsistent direction parameter usage")


if __name__ == "__main__":
    main()
