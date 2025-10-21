"""Step 1: Comprehensive diagnostic of multiple scattering integration.

This script checks:
1. Integration grid (U, V) values
2. Radiation mask coverage
3. Integrand magnitudes (Ikc, Ic)
4. Propagator values (Fp, Fm, Gp, Gm)
5. Kirchhoff-complementary terms (K1, K2, K3)
6. Complementary terms (gc1-gc14)

Results are categorized by:
- Surface roughness ratio (ℓ/σ)
- Single scattering vs multiple scattering
- Polarization (VV, HH, HV)
"""

import numpy as np
from pathlib import Path
from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState
from mwrtms.scattering.surface.iem.multiple_scattering import (
    _prepare_geometry_params,
    _build_quadrature,
    PhysicsParams,
    SurfaceParams,
    _make_Wn_provider,
    _assemble_integrands,
    _build_propagators,
    _build_gkc1,
    _build_gkc2,
    _build_gkc3,
    _build_gc_block1,
    _build_gc_block2,
)

# Load NMM3D data
lut_path = Path("data/NMM3D_LUT_NRCS_40degree.dat")
table = np.loadtxt(lut_path)
mask = np.isclose(table[:, 0], 40.0)
rows = table[mask]

frequency_ghz = 5.405
wavelength_m = 0.3 / frequency_ghz
k = 2 * np.pi / wavelength_m
radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=40.0)

# Select test cases: one per ratio with finite HV and moderate roughness
test_cases = []
for target_ratio in [4.0, 7.0, 10.0, 15.0]:
    for row in rows:
        theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
        ks = k * rms_norm * wavelength_m
        if np.isclose(ratio, target_ratio) and np.isfinite(hv_ref) and 0.5 < ks < 1.0:
            test_cases.append(row)
            break

print("=" * 100)
print("STEP 1: MULTIPLE SCATTERING INTEGRATION DIAGNOSTICS")
print("=" * 100)

for case_idx, row in enumerate(test_cases):
    theta, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
    
    sigma = rms_norm * wavelength_m
    corr_len = ratio * sigma
    rms_height_cm = sigma * 100.0
    correlation_length_cm = corr_len * 100.0
    soil_permittivity = complex(float(eps_r), float(eps_i))
    
    ks = k * sigma
    kl = k * corr_len
    
    print(f"\n{'=' * 100}")
    print(f"TEST CASE {case_idx + 1}: ℓ/σ = {ratio:.1f}")
    print(f"{'=' * 100}")
    print(f"Parameters: εr={eps_r:.1f}-j{eps_i:.1f}, kσ={ks:.3f}, kℓ={kl:.3f}")
    print(f"NMM3D: VV={vv_ref:.2f} dB, HH={hh_ref:.2f} dB, HV={hv_ref:.2f} dB")
    
    # ========================================================================
    # SINGLE SCATTERING
    # ========================================================================
    print(f"\n{'-' * 100}")
    print("SINGLE SCATTERING (Current AIEM Implementation)")
    print(f"{'-' * 100}")
    
    vv_single = mwRTMs.compute_soil_backscatter(
        model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity, correlation='exponential',
        polarization=PolarizationState.VV,
    )
    hh_single = mwRTMs.compute_soil_backscatter(
        model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity, correlation='exponential',
        polarization=PolarizationState.HH,
    )
    hv_single = mwRTMs.compute_soil_backscatter(
        model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity, correlation='exponential',
        polarization=PolarizationState.HV,
    )
    
    vv_single_db = 10*np.log10(vv_single) if vv_single > 0 else -999
    hh_single_db = 10*np.log10(hh_single) if hh_single > 0 else -999
    hv_single_db = 10*np.log10(hv_single) if hv_single > 0 else -999
    
    print(f"VV: {vv_single_db:7.2f} dB (σ⁰={vv_single:.3e}) | Error: {vv_single_db-vv_ref:+6.2f} dB")
    print(f"HH: {hh_single_db:7.2f} dB (σ⁰={hh_single:.3e}) | Error: {hh_single_db-hh_ref:+6.2f} dB")
    print(f"HV: {hv_single_db:7.2f} dB (σ⁰={hv_single:.3e}) | Error: {hv_single_db-hv_ref:+6.2f} dB")
    
    # ========================================================================
    # MULTIPLE SCATTERING DIAGNOSTICS
    # ========================================================================
    print(f"\n{'-' * 100}")
    print("MULTIPLE SCATTERING DIAGNOSTICS")
    print(f"{'-' * 100}")
    
    # Setup geometry and parameters
    theta_i = np.radians(40.0)
    theta_s = theta_i  # Backscatter
    phi_i = 0.0
    phi_s = np.pi  # Backscatter convention
    
    geom = _prepare_geometry_params(theta_i, theta_s, phi_i, phi_s, k)
    phys = PhysicsParams(k=k, er=soil_permittivity)
    surf = SurfaceParams(type='exponential', ks=ks, kl=kl, sigma=sigma, corr_length_m=corr_len)
    
    # Test with different grid sizes
    for n_points in [33, 65]:
        print(f"\n  Grid Size: {n_points} × {n_points} points, nmax=6")
        print(f"  {'-' * 96}")
        
        quad = _build_quadrature(surf, k, n_points=n_points, nmax=6)
        
        # Integration domain
        U = quad.U
        V = quad.V
        umax = 5.0 / max(surf.corr_length_m, 1e-6)
        print(f"  Integration domain: U,V ∈ [{-umax:.3f}, {+umax:.3f}]")
        print(f"  Grid spacing: ΔU = ΔV = {U[0,1]-U[0,0]:.4f}")
        
        # Compute wavenumbers
        q1 = np.sqrt(np.maximum(k**2 - (U**2 + V**2), 0.0))
        q2 = np.sqrt(phys.er * k**2 - (U**2 + V**2))
        
        print(f"  q1 (air):      min={q1.min():.2f}, max={q1.max():.2f}, mean={q1.mean():.2f}")
        print(f"  |q2| (soil):   min={np.abs(q2).min():.2f}, max={np.abs(q2).max():.2f}")
        
        # Radiation condition
        qmin = 1e-6
        rad = (np.real(q1) > qmin) | (np.real(q2) > qmin)
        print(f"  Radiation mask: {rad.sum()}/{rad.size} points valid ({100*rad.sum()/rad.size:.1f}%)")
        
        # Roughness spectrum
        wn_provider = _make_Wn_provider(surf)
        Wn1 = wn_provider(U, V, 1)
        print(f"  W₁(u,v): min={Wn1.min():.3e}, max={Wn1.max():.3e}, mean={Wn1.mean():.3e}")
        
        # Test HV polarization
        print(f"\n  HV Polarization Analysis:")
        print(f"  {'-' * 96}")
        
        # Build propagators
        propagators = _build_propagators(U, V, q1, q2, k, phys.er, geom, 'hv')
        
        for prop_name in ['Fp', 'Fm', 'Gp', 'Gm']:
            P = propagators[prop_name]
            print(f"    {prop_name}: |min|={np.abs(P).min():.3e}, |max|={np.abs(P).max():.3e}, "
                  f"|mean|={np.abs(P).mean():.3e}")
        
        # Build Kirchhoff-complementary terms
        K1 = _build_gkc1(U, V, geom, q1, surf, wn_provider, quad.Nmax)
        K2 = _build_gkc2(U, V, geom, q1, surf, wn_provider, quad.Nmax)
        K3 = _build_gkc3(U, V, geom, q1, surf, wn_provider, quad.Nmax)
        
        print(f"    K1: |min|={np.abs(K1).min():.3e}, |max|={np.abs(K1).max():.3e}, "
              f"|mean|={np.abs(K1).mean():.3e}")
        print(f"    K2: |min|={np.abs(K2).min():.3e}, |max|={np.abs(K2).max():.3e}, "
              f"|mean|={np.abs(K2).mean():.3e}")
        print(f"    K3: |min|={np.abs(K3).min():.3e}, |max|={np.abs(K3).max():.3e}, "
              f"|mean|={np.abs(K3).mean():.3e}")
        
        # Build complementary terms
        C1 = _build_gc_block1(U, V, geom, q1, q1, surf, wn_provider, quad.Nmax)
        C2 = _build_gc_block2(U, V, geom, q1, q1, surf, wn_provider, quad.Nmax)
        
        # Sample a few gc terms
        for gc_name in ['gc1', 'gc5', 'gc9']:
            if gc_name in C1:
                gc = C1[gc_name]
            else:
                gc = C2[gc_name]
            print(f"    {gc_name}: |min|={np.abs(gc).min():.3e}, |max|={np.abs(gc).max():.3e}, "
                  f"|mean|={np.abs(gc).mean():.3e}")
        
        # Compute integrands
        integrand_kc, integrand_c = _assemble_integrands(
            U, V, q1, q2, k, phys.er, geom, surf, wn_provider, quad.Nmax, 'hv'
        )
        
        print(f"    Integrand KC: |min|={np.abs(integrand_kc).min():.3e}, "
              f"|max|={np.abs(integrand_kc).max():.3e}, |mean|={np.abs(integrand_kc).mean():.3e}")
        print(f"    Integrand C:  |min|={np.abs(integrand_c).min():.3e}, "
              f"|max|={np.abs(integrand_c).max():.3e}, |mean|={np.abs(integrand_c).mean():.3e}")
        
        # Apply radiation mask and compute integrals
        Ikc = np.real(integrand_kc) * rad
        Ic = np.real(integrand_c) * rad
        
        print(f"    After rad mask - KC: non-zero={np.count_nonzero(Ikc)}/{Ikc.size}, "
              f"sum={np.sum(Ikc):.3e}")
        print(f"    After rad mask - C:  non-zero={np.count_nonzero(Ic)}/{Ic.size}, "
              f"sum={np.sum(Ic):.3e}")
        
        # Compute weights and final integral
        W2D = np.outer(quad.wu, quad.wv)
        val_kc = (k**2 / (8.0 * np.pi)) * np.sum(Ikc * W2D)
        val_c = (k**2 / (64.0 * np.pi)) * np.sum(Ic * W2D)
        val_total = val_kc + val_c
        
        print(f"    KC contribution:    {val_kc:.6e}")
        print(f"    C contribution:     {val_c:.6e}")
        print(f"    Total MS:           {val_total:.6e}")
        
        if val_total > 0:
            val_total_db = 10*np.log10(val_total)
            print(f"    Total MS (dB):      {val_total_db:.2f} dB")
        else:
            print(f"    Total MS (dB):      -inf (zero or negative)")
    
    # ========================================================================
    # FULL COMPUTATION WITH MS
    # ========================================================================
    print(f"\n{'-' * 100}")
    print("FULL AIEM WITH MULTIPLE SCATTERING")
    print(f"{'-' * 100}")
    
    hv_with_ms = mwRTMs.compute_soil_backscatter(
        model='aiem', radar_config=radar_config, frequency_ghz=frequency_ghz,
        rms_height_cm=rms_height_cm, correlation_length_cm=correlation_length_cm,
        soil_permittivity=soil_permittivity, correlation='exponential',
        polarization=PolarizationState.HV,
        include_multiple_scattering=True,
        ms_quadrature_points=65,
        ms_spectral_terms=6,
    )
    
    hv_with_ms_db = 10*np.log10(hv_with_ms) if hv_with_ms > 0 else -999
    ms_contrib = hv_with_ms - hv_single
    
    print(f"HV Single:      {hv_single_db:7.2f} dB (σ⁰={hv_single:.3e})")
    print(f"HV with MS:     {hv_with_ms_db:7.2f} dB (σ⁰={hv_with_ms:.3e})")
    print(f"MS Contribution: {ms_contrib:.3e} (Δ={hv_with_ms_db-hv_single_db:.2f} dB)")
    print(f"NMM3D HV:       {hv_ref:7.2f} dB")
    print(f"Error:          {hv_with_ms_db-hv_ref:+6.2f} dB")

print(f"\n{'=' * 100}")
print("DIAGNOSTIC COMPLETE")
print(f"{'=' * 100}")
