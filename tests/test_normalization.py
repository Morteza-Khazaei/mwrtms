#!/usr/bin/env python3
"""Test different normalization factors."""

import numpy as np

# Test the magnitude of (2π)^n for different n
print("Normalization factor magnitudes:")
print("=" * 50)

for n in [1.5, 2, 3, 4, 5, 10]:
    factor = (2.0 * np.pi) ** n
    print(f"(2π)^{n:4.1f} = {factor:.6e}")

print()
print("Analysis:")
print("=" * 50)
print("(2π)^10 ≈ 9.35 × 10^6 is HUGE!")
print()
print("In the spectrum W(u,v,n), this gets:")
print("- Raised to power n in series (up to n=8)")
print("- Multiplied across two series sums")
print("- Integrated over 129×129 grid")
print()
print("This could cause:")
print("1. Numerical overflow in intermediate calculations")
print("2. Loss of precision in floating point")
print("3. Integration producing inf or nan")
print()
print("The spectrum value at each point:")
factor_10 = (2.0 * np.pi) ** 10
sigma2 = (0.005) ** 2  # 0.5 cm in meters
kl = 2.264
n = 1
W_magnitude = factor_10 * sigma2 * (kl / n) ** 2
print(f"W(u,v,1) ~ {factor_10:.2e} × {sigma2:.2e} × {(kl/n)**2:.2f}")
print(f"        ~ {W_magnitude:.2e}")
print()
print("When raised to power 8 in series:")
print(f"W^8 ~ ({W_magnitude:.2e})^8 = {W_magnitude**8:.2e}")
print()
print("This is astronomical and will cause numerical issues!")
