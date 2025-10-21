"""Test different normalization factors for exponential spectrum."""

import numpy as np

# Standard 2D Fourier transform of exp(-r/L) is:
# W(κx, κy) = ∫∫ exp(-r/L) exp(-i(κx×x + κy×y)) dx dy
# where r = sqrt(x² + y²), κ = sqrt(κx² + κy²)
#
# Result in polar coordinates:
# W(κ) = 2πL² / (1 + (Lκ)²)^(3/2)
#
# For nth power: exp(-nr/L) → equivalent to L→L/n
# W^(n)(κ) = 2π(L/n)² / (1 + (Lκ/n)²)^(3/2)

print("Standard 2D Fourier Transform Results:")
print("=" * 70)

# For reference
L = 0.06  # m
n_vals = [1, 2, 3, 4]

print(f"\nFor exp(-r/L) with L = {L*100:.1f} cm:\n")

for n in n_vals:
    # Standard FT result
    prefactor_standard = 2 * np.pi * (L/n)**2

    # Alternative: might need to account for σ² separately
    # or different conventions

    print(f"n = {n}:")
    print(f"  Standard FT prefactor: 2π(L/n)² = {prefactor_standard:.6e}")
    print(f"  Which equals: 2πL²/n² = {2*np.pi*L**2/n**2:.6e}")

    # Can also write as:
    # W^(n)(κ) = 2πL²/n² / [1 + (Lκ/n)²]^(3/2)
    #          = 2πL² / [n²(1 + (Lκ/n)²)]^(3/2)
    #          = 2πL² / [n² + (Lκ)²]^(3/2)  ← if we multiply num/denom by n³

    # Actually, let's be more careful:
    # [n²(1 + (Lκ/n)²)]^(3/2) = [n² + (Lκ)²]^(3/2) / n³ × n³
    # NO! Let me recompute:
    # n²(1 + (Lκ/n)²) = n² + (Lκ)²
    # So [n²(1 + (Lκ/n)²)]^(3/2) = [n² + (Lκ)²]^(3/2)

    # Therefore:
    # W^(n)(κ) = 2πL²/n² / [1 + (Lκ/n)²]^(3/2)
    #          = 2πL² / n² × [n² + (Lκ)²]^(-3/2) × n³
    #          = 2πL²n / [n² + (Lκ)²]^(3/2)

    prefactor_rewritten = 2 * np.pi * L**2 * n
    print(f"  Rewritten: 2πL²n / [n² + (Lκ)²]^(3/2) prefactor = {prefactor_rewritten:.6e}")
    print()

print("\n" + "=" * 70)
print("CONCLUSION:")
print("=" * 70)
print("\nThe correct formula should be:")
print("  W^(n)(κ) = 2πL²n / [n² + (Lκ)²]^(3/2) × σ²")
print("\nNOT:")
print("  W^(n)(κ) = 4πL²n / [n² + (Lκ)²]^(3/2) × σ²")
print("\nThe factor should be 2π, not 4π!")
print("\nThis would give a 2x reduction, not the 10^7 needed...")
print("\nSo this alone doesn't explain the discrepancy.")
