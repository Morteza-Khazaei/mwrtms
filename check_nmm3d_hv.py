"""Check NMM3D HV values to find finite cases."""

import numpy as np
from pathlib import Path

# Load NMM3D data
lut_path = Path("data/NMM3D_LUT_NRCS_40degree.dat")
table = np.loadtxt(lut_path)
mask = np.isclose(table[:, 0], 40.0)
rows = table[mask]

print("=" * 80)
print("NMM3D HV Value Analysis")
print("=" * 80)

# Check HV column (index 7)
hv_values = rows[:, 7]

print(f"\nTotal samples: {len(hv_values)}")
print(f"Finite HV values: {np.isfinite(hv_values).sum()}")
print(f"Infinite HV values: {np.isinf(hv_values).sum()}")
print(f"NaN HV values: {np.isnan(hv_values).sum()}")

finite_mask = np.isfinite(hv_values)
if finite_mask.any():
    finite_hv = hv_values[finite_mask]
    print(f"\nFinite HV statistics:")
    print(f"  Min: {finite_hv.min():.2f} dB")
    print(f"  Max: {finite_hv.max():.2f} dB")
    print(f"  Mean: {finite_hv.mean():.2f} dB")
    print(f"  Median: {np.median(finite_hv):.2f} dB")
    
    print(f"\nFirst 10 samples with finite HV:")
    print(f"{'Idx':<5} {'ℓ/σ':<6} {'εr':<6} {'kσ':<8} {'VV (dB)':<10} {'HH (dB)':<10} {'HV (dB)':<10}")
    print("-" * 80)
    
    frequency_ghz = 5.405
    wavelength_m = 0.3 / frequency_ghz
    k = 2 * np.pi / wavelength_m
    
    count = 0
    for i, row in enumerate(rows):
        if np.isfinite(row[7]) and count < 10:
            theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
            ks = k * rms_norm * wavelength_m
            print(f"{i:<5} {ratio:<6.1f} {eps_r:<6.1f} {ks:<8.3f} {vv_ref:<10.2f} {hh_ref:<10.2f} {hv_ref:<10.2f}")
            count += 1
else:
    print("\nNO FINITE HV VALUES FOUND!")
    print("This means NMM3D dataset has -inf for all HV values.")
    print("The test comparison is meaningless for HV polarization.")

print("\n" + "=" * 80)
