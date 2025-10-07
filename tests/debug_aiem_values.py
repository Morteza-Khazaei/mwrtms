"""Debug script to check AIEM intermediate values."""

import numpy as np
from src.mwrtms.scattering.iem.fresnel_utils import compute_fresnel_incident, compute_fresnel_specular
from src.mwrtms.scattering.iem.transition import compute_transition_function
from src.mwrtms.scattering.iem.spectrum_aiem import compute_aiem_spectrum
from src.mwrtms.scattering.iem.geometry_utils import compute_spatial_frequency

# Test case from NMM3D: ratio=4, eps_r=3+1j, rms_norm=0.021
frequency_ghz = 5.405
lam = 0.3 / frequency_ghz  # wavelength in meters
k = 2.0 * np.pi / lam

theta_i = np.deg2rad(40.0)
theta_s = theta_i  # monostatic
phi_s = np.pi  # backscatter

eps_r = 3.0 + 1.0j
ratio = 4.0
rms_norm = 0.021

sigma = rms_norm * lam
corr_len = ratio * sigma
ks = k * sigma
kl = k * corr_len

print("="*60)
print("AIEM Intermediate Values Debug")
print("="*60)
print(f"\nInput parameters:")
print(f"  Frequency: {frequency_ghz} GHz")
print(f"  Wavelength: {lam*100:.4f} cm")
print(f"  k: {k:.4f} rad/m")
print(f"  θ_i: {np.rad2deg(theta_i):.1f}°")
print(f"  ε_r: {eps_r}")
print(f"  ℓ/σ: {ratio}")
print(f"  rms_norm: {rms_norm}")
print(f"  σ: {sigma*100:.4f} cm")
print(f"  ℓ: {corr_len*100:.4f} cm")
print(f"  ks: {ks:.4f}")
print(f"  kl: {kl:.4f}")

# Compute Fresnel coefficients
print(f"\n{'='*60}")
print("Fresnel Coefficients")
print("="*60)

Rvi, Rhi, Rvhi = compute_fresnel_incident(eps_r, theta_i)
print(f"At incident angle:")
print(f"  Rv: {Rvi:.6f}")
print(f"  Rh: {Rhi:.6f}")
print(f"  |Rv|: {np.abs(Rvi):.6f}")
print(f"  |Rh|: {np.abs(Rhi):.6f}")

Rvl, Rhl, Rvhl = compute_fresnel_specular(eps_r, theta_i, theta_s, phi_s)
print(f"\nAt specular angle:")
print(f"  Rv: {Rvl:.6f}")
print(f"  Rh: {Rhl:.6f}")
print(f"  |Rv|: {np.abs(Rvl):.6f}")
print(f"  |Rh|: {np.abs(Rhl):.6f}")

# Compute spectrum
print(f"\n{'='*60}")
print("Roughness Spectrum")
print("="*60)

K = compute_spatial_frequency(kl, theta_i, theta_s, phi_s, 0.0)
print(f"K: {K:.4f}")

print(f"\nW^(n)(K) for n=1..5:")
for n in range(1, 6):
    W_n = compute_aiem_spectrum(kl, K, n, 'exponential')
    print(f"  n={n}: W^({n}) = {W_n:.6e}")

# Compute transition function
print(f"\n{'='*60}")
print("Transition Function")
print("="*60)

cs = np.cos(theta_i)
n_terms = 10
spectra = np.array([compute_aiem_spectrum(kl, K, n, 'exponential') for n in range(1, n_terms+1)])

Tfv, Tfh = compute_transition_function(eps_r, theta_i, ks, cs, spectra, n_terms)
print(f"Tfv: {Tfv:.6f}")
print(f"Tfh: {Tfh:.6f}")

# Transitioned Fresnel coefficients
Rvtran = Rvi + (Rvl - Rvi) * Tfv
Rhtran = Rhi + (Rhl - Rhi) * Tfh
print(f"\nTransitioned coefficients:")
print(f"  Rv_tran: {Rvtran:.6f}")
print(f"  Rh_tran: {Rhtran:.6f}")
print(f"  |Rv_tran|: {np.abs(Rvtran):.6f}")
print(f"  |Rh_tran|: {np.abs(Rhtran):.6f}")

print(f"\n{'='*60}")
print("Expected NMM3D values (from LUT):")
print("="*60)
print(f"  VV: -27.29 dB")
print(f"  HH: -28.25 dB")
print(f"  HV: -Inf dB")
