"""Diagnose HV bias dependence on roughness ratio."""

import numpy as np
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry, PolarizationState
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Test parameters
freq = 5.3e9  # Hz
theta_deg = 40.0
er = 12.0 + 1.8j

# Test different roughness ratios
test_cases = [
    # (sigma, L, ratio)
    (0.010, 0.040, 4.0),   # Rough: ℓ/σ = 4
    (0.015, 0.060, 4.0),   # Rough: ℓ/σ = 4
    (0.010, 0.070, 7.0),   # Medium: ℓ/σ = 7
    (0.015, 0.105, 7.0),   # Medium: ℓ/σ = 7
    (0.010, 0.100, 10.0),  # Smooth: ℓ/σ = 10
    (0.015, 0.150, 10.0),  # Smooth: ℓ/σ = 10
    (0.010, 0.150, 15.0),  # Very smooth: ℓ/σ = 15
    (0.015, 0.225, 15.0),  # Very smooth: ℓ/σ = 15
]

print("\n" + "="*90)
print("HV CROSS-POL ROUGHNESS DEPENDENCE ANALYSIS")
print("="*90)

results = []

wave = ElectromagneticWave(freq)
geometry = ScatteringGeometry(theta_i_deg=theta_deg)
air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(er)

for sigma, L, ratio in test_cases:
    surface = build_surface_from_statistics(sigma, L, correlation_type="exponential")

    # Model with MS
    model_ms = AIEMModel(
        wave, geometry, surface,
        correlation_type="exponential",
        include_multiple_scattering=True,
        enable_guardrails=False,
        ms_quadrature_points=129,
        ms_spectral_terms=8,
    )

    # Model without MS
    model_no_ms = AIEMModel(
        wave, geometry, surface,
        correlation_type="exponential",
        include_multiple_scattering=False,
        enable_guardrails=False,
    )

    # Compute
    hv_total = model_ms.compute(air, soil, PolarizationState.HV)
    hv_single = model_no_ms.compute(air, soil, PolarizationState.HV)
    hv_ms = hv_total - hv_single

    vv_total = model_ms.compute(air, soil, PolarizationState.VV)
    hh_total = model_ms.compute(air, soil, PolarizationState.HH)

    ks = wave.wavenumber * sigma
    kl = wave.wavenumber * L

    results.append({
        'ratio': ratio,
        'sigma_cm': sigma * 100,
        'L_cm': L * 100,
        'ks': ks,
        'kl': kl,
        'hv_single': hv_single,
        'hv_ms': hv_ms,
        'hv_total': hv_total,
        'hv_total_db': 10*np.log10(hv_total) if hv_total > 0 else -999,
        'hv_ms_db': 10*np.log10(hv_ms) if hv_ms > 0 else -999,
        'ms_fraction': hv_ms / hv_total if hv_total > 0 else 0,
        'vv_total': vv_total,
        'hh_total': hh_total,
        'hv_vv_ratio': hv_total / vv_total if vv_total > 0 else 0,
        'hv_hh_ratio': hv_total / hh_total if hh_total > 0 else 0,
    })

print(f"\nRoughness Ratio Analysis:")
print(f"{'Ratio':<8} {'σ(cm)':<8} {'L(cm)':<8} {'kσ':<8} {'kL':<8} {'HV(dB)':<10} {'HV_MS(dB)':<12} {'MS%':<8} {'HV/VV':<8}")
print("-" * 90)

for row in results:
    print(f"{row['ratio']:<8.1f} {row['sigma_cm']:<8.2f} {row['L_cm']:<8.2f} "
          f"{row['ks']:<8.3f} {row['kl']:<8.3f} {row['hv_total_db']:<10.2f} "
          f"{row['hv_ms_db']:<12.2f} {row['ms_fraction']*100:<8.1f} "
          f"{10*np.log10(row['hv_vv_ratio']):<8.2f}")

# Group by ratio and compute statistics
print(f"\n" + "="*90)
print("SUMMARY BY ROUGHNESS RATIO:")
print("="*90)

# Collect by ratio
ratio_groups = {}
for row in results:
    r = row['ratio']
    if r not in ratio_groups:
        ratio_groups[r] = []
    ratio_groups[r].append(row)

for ratio in sorted(ratio_groups.keys()):
    subset = ratio_groups[ratio]
    hv_dbs = [r['hv_total_db'] for r in subset]
    ms_dbs = [r['hv_ms_db'] for r in subset]
    ms_fracs = [r['ms_fraction'] for r in subset]
    hv_vv_ratios = [10*np.log10(r['hv_vv_ratio']) for r in subset]
    kss = [r['ks'] for r in subset]
    kls = [r['kl'] for r in subset]

    print(f"\nℓ/σ = {ratio:.1f}:")
    print(f"  HV total (dB):     {np.mean(hv_dbs):.2f} ± {np.std(hv_dbs):.2f}")
    print(f"  HV from MS (dB):   {np.mean(ms_dbs):.2f} ± {np.std(ms_dbs):.2f}")
    print(f"  MS contributes:    {np.mean(ms_fracs)*100:.1f}% of total")
    print(f"  HV/VV (dB):        {np.mean(hv_vv_ratios):.2f}")
    print(f"  Mean kσ:           {np.mean(kss):.3f}")
    print(f"  Mean kL:           {np.mean(kls):.3f}")

print(f"\n" + "="*90)
print("\nKEY OBSERVATIONS:")
print("-" * 90)

# Check for roughness-dependent trends
hv_4 = [r['hv_total_db'] for r in results if r['ratio'] == 4.0]
hv_15 = [r['hv_total_db'] for r in results if r['ratio'] == 15.0]
ratio_4 = np.mean(hv_4)
ratio_15 = np.mean(hv_15)
trend = ratio_4 - ratio_15

print(f"1. HV bias vs roughness trend: {trend:.2f} dB (rough→smooth)")
if abs(trend) > 10:
    print(f"   ⚠️  STRONG roughness dependence detected!")
elif abs(trend) > 5:
    print(f"   ⚠️  Moderate roughness dependence")
else:
    print(f"   ✓ Weak roughness dependence")

# Check MS contribution
ms_fracs_4 = [r['ms_fraction'] for r in results if r['ratio'] == 4.0]
ms_fracs_15 = [r['ms_fraction'] for r in results if r['ratio'] == 15.0]
ms_contrib_4 = np.mean(ms_fracs_4)
ms_contrib_15 = np.mean(ms_fracs_15)
print(f"\n2. MS contribution:")
print(f"   Rough (ℓ/σ=4):   {ms_contrib_4*100:.1f}%")
print(f"   Smooth (ℓ/σ=15): {ms_contrib_15*100:.1f}%")

# Check if abs() is the issue
print(f"\n3. Suggested diagnosis:")
if ratio_4 > -5:  # Too large for rough surfaces
    print(f"   - Rough surfaces (ℓ/σ=4): HV too large → likely over-integration")
    print(f"   - Possible cause: np.abs() preventing destructive interference")
if ratio_15 < -10:  # Too small for smooth surfaces
    print(f"   - Smooth surfaces (ℓ/σ=15): HV too small → likely under-integration")
    print(f"   - Possible cause: Integration domain too narrow or excessive cancellation")
