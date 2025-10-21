"""Debug spectrum magnitude with and without (2π)^10 factor."""

import numpy as np

# Test parameters (typical)
freq = 5.3e9  # Hz
c = 3e8
wavelength = c / freq
k = 2 * np.pi / wavelength

sigma = 0.015  # 1.5 cm
L = 0.06  # 6 cm
ks = k * sigma
kl = k * L
sigma2 = sigma**2

print(f"\nSpectrum Magnitude Comparison")
print(f"=" * 70)
print(f"Parameters: σ={sigma*100:.1f} cm, L={L*100:.1f} cm")
print(f"  k = {k:.2f} m⁻¹")
print(f"  kσ = {ks:.3f}")
print(f"  kL = {kl:.3f}")

# Test at typical spectral point
u = 2.0 * k  # Typical integration point
v = 0.0
kappa = np.sqrt(u**2 + v**2)

print(f"\nAt typical point (U={u/k:.1f}k, V=0):")
print(f"  κ = {kappa:.2f} m⁻¹")

# Test for different orders n
for n in [1, 2, 3, 4]:
    # OLD formula with (2π)^10
    denom_old = 1.0 + ((kl * kappa) / n) ** 2
    W_old = (2.0 * np.pi)**10 * sigma2 * (kl / n) ** 2 * denom_old ** (-1.5)

    # NEW formula (correct)
    numerator_new = 4.0 * np.pi * n * kl**2
    denominator_new = (n**2 + (kl * kappa)**2) ** 1.5
    W_new = sigma2 * numerator_new / denominator_new

    ratio = W_old / W_new if W_new > 0 else float('inf')

    print(f"\nn = {n}:")
    print(f"  W^({n}) OLD (with (2π)^10): {W_old:.6e}")
    print(f"  W^({n}) NEW (correct):       {W_new:.6e}")
    print(f"  Ratio (OLD/NEW):            {ratio:.2e}")

# Check what the (2π)^10 factor is
factor = (2.0 * np.pi)**10
print(f"\n" + "=" * 70)
print(f"(2π)^10 = {factor:.6e}")
print(f"\nThis factor is scaling the spectrum by ~{factor:.0f}x!")
