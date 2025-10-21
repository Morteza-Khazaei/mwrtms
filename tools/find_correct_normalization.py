"""Find the correct normalization factor for exponential spectrum."""

import numpy as np

# The ratio we found
ratio_found = 7.63e6

# (2π)^10
factor_old = (2.0 * np.pi)**10
print(f"(2π)^10 = {factor_old:.6e}")
print(f"Ratio OLD/NEW at test point = {ratio_found:.6e}")

# What factor would we need?
print(f"\nIf we need to scale by ~{ratio_found:.2e}, what could this be?")

# Test various powers of 2π
for p in range(1, 15):
    val = (2.0 * np.pi)**p
    if abs(val - ratio_found) / ratio_found < 0.5:
        print(f"  (2π)^{p} = {val:.6e} ← Close match!")

# Test combinations
print(f"\nOther factors to test:")
print(f"  4π = {4*np.pi:.6e}")
print(f"  (4π)² = {(4*np.pi)**2:.6e}")
print(f"  (4π)³ = {(4*np.pi)**3:.6e}")
print(f"  (4π)⁴ = {(4*np.pi)**4:.6e}")
print(f"  (4π)⁵ = {(4*np.pi)**5:.6e}")

# Check if it's exactly (2π)^10
print(f"\n(2π)^10 / ratio = {factor_old / ratio_found:.3f}")

# Check sqrt or other factors
print(f"\nCheck specific values:")
print(f"  (2π)^10 / (4π) = {factor_old / (4*np.pi):.6e}")
print(f"  (2π)^8 * 4π = {(2*np.pi)**8 * 4 * np.pi:.6e}")
print(f"  (4π)^5 = {(4*np.pi)**5:.6e}")

# The empirical factor was (2π)^10
# The NEW formula has 4π in numerator
# So effective scaling relative to NEW is (2π)^10 / (4π) ?
effective_with_4pi = factor_old / (4.0 * np.pi)
print(f"\n(2π)^10 / (4π) = {effective_with_4pi:.6e}")
print(f"This is (2π)^10 / (4π) = {effective_with_4pi:.6e}")

# Actually the formula has:
# NEW: 4πnℓ² / [...]
# OLD: (2π)^10 * (ℓ/n)² / [...]
# At n=1: OLD has (2π)^10 * ℓ²
# At n=1: NEW has 4π * 1 * ℓ² = 4πℓ²
# Ratio: (2π)^10 / (4π) = (2π)^10 / (4π) = (2π)^9 / 2

print(f"\n(2π)^9 / 2 = {(2*np.pi)**9 / 2:.6e}")
print(f"(2π)^9 * 2 = {(2*np.pi)**9 * 2:.6e}")
