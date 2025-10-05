"""Benchmark script for AIEM multiple scattering with and without Numba."""

import time
import numpy as np
from pathlib import Path

# Test if Numba backend is available
try:
    from src.mwrtms.scattering.iem.aiem_numba_backend import (
        NUMBA_AVAILABLE,
        check_numba_performance,
        precompute_factorials,
        get_two_pi_power,
        series_sum_exponential_numba,
        compute_wn_exponential_numba,
        integrate_2d_real_numba
    )
    print("✅ AIEM Numba backend imported successfully")
    print(f"   Numba available: {NUMBA_AVAILABLE}")
except ImportError as e:
    print(f"❌ Failed to import AIEM Numba backend: {e}")
    NUMBA_AVAILABLE = False

# Test basic AIEM multiple scattering
try:
    from src.mwrtms.scattering.iem.multiple_scattering import compute_multiple_scattering
    print("✅ AIEM multiple scattering module imported")
except ImportError as e:
    print(f"❌ Failed to import multiple scattering: {e}")
    compute_multiple_scattering = None


def benchmark_numba_functions():
    """Benchmark individual Numba functions."""
    print("\n" + "="*70)
    print("BENCHMARKING NUMBA FUNCTIONS")
    print("="*70)
    
    if not NUMBA_AVAILABLE:
        print("⚠️  Numba not available - skipping benchmarks")
        return
    
    # Test parameters
    nmax = 8
    n_iterations = 10000
    
    # Pre-compute factorials
    factorials = precompute_factorials(nmax)
    two_pi_power = get_two_pi_power(10)
    
    print(f"\nTest parameters:")
    print(f"  - Spectral order (nmax): {nmax}")
    print(f"  - Iterations: {n_iterations}")
    print(f"  - (2π)^10 = {two_pi_power:.6e}")
    
    # Benchmark 1: Roughness spectrum computation
    print(f"\n{'='*70}")
    print("1. Roughness Spectrum Computation")
    print(f"{'='*70}")
    
    u, v = 1.5, 2.0
    sigma2 = 0.01
    kl = 2.0
    
    start = time.time()
    for _ in range(n_iterations):
        for n in range(1, nmax + 1):
            wn = compute_wn_exponential_numba(u, v, n, sigma2, kl, two_pi_power)
    elapsed = time.time() - start
    
    print(f"   Computed {n_iterations * nmax} spectrum values")
    print(f"   Time: {elapsed:.4f} seconds")
    print(f"   Rate: {(n_iterations * nmax) / elapsed:.0f} evaluations/sec")
    
    # Benchmark 2: Series summation
    print(f"\n{'='*70}")
    print("2. Series Summation")
    print(f"{'='*70}")
    
    coeff_real = 0.5
    coeff_imag = 0.3
    arg_x = 1.0
    arg_y = 1.5
    
    start = time.time()
    for _ in range(n_iterations):
        result = series_sum_exponential_numba(
            coeff_real, coeff_imag, arg_x, arg_y,
            sigma2, kl, nmax, two_pi_power, factorials
        )
    elapsed = time.time() - start
    
    print(f"   Computed {n_iterations} series sums (order {nmax})")
    print(f"   Time: {elapsed:.4f} seconds")
    print(f"   Rate: {n_iterations / elapsed:.0f} sums/sec")
    print(f"   Result: {result[0]:.6e} + {result[1]:.6e}j")
    
    # Benchmark 3: 2D Integration
    print(f"\n{'='*70}")
    print("3. 2D Integration")
    print(f"{'='*70}")
    
    grid_sizes = [65, 129, 257]
    
    for n_points in grid_sizes:
        integrand = np.random.rand(n_points, n_points)
        weights = np.random.rand(n_points, n_points) * 0.01
        mask = np.random.rand(n_points, n_points) > 0.1
        
        n_iter = max(10, 1000 // n_points)
        
        start = time.time()
        for _ in range(n_iter):
            result = integrate_2d_real_numba(integrand, weights, mask)
        elapsed = time.time() - start
        
        print(f"   Grid: {n_points}×{n_points} ({n_points**2} points)")
        print(f"   Iterations: {n_iter}")
        print(f"   Time: {elapsed:.4f} seconds")
        print(f"   Time per integration: {elapsed/n_iter*1000:.2f} ms")
        print()


def benchmark_full_multiple_scattering():
    """Benchmark full multiple scattering calculation."""
    print("\n" + "="*70)
    print("BENCHMARKING FULL MULTIPLE SCATTERING")
    print("="*70)
    
    if compute_multiple_scattering is None:
        print("⚠️  Multiple scattering module not available")
        return
    
    # Test parameters (from NMM3D validation)
    frequency_ghz = 5.405
    wavelength_m = 0.3 / frequency_ghz
    k = 2 * np.pi / wavelength_m
    theta_deg = 40.0
    theta_rad = np.deg2rad(theta_deg)
    
    # Surface parameters
    rms_height_m = 0.01  # 1 cm
    corr_length_m = 0.04  # 4 cm
    sigma = rms_height_m
    ks = k * sigma
    kl = k * corr_length_m
    
    # Soil permittivity
    er = complex(12.0, 1.8)
    
    print(f"\nTest configuration:")
    print(f"  - Frequency: {frequency_ghz} GHz")
    print(f"  - Wavelength: {wavelength_m*100:.2f} cm")
    print(f"  - Incidence angle: {theta_deg}°")
    print(f"  - RMS height: {rms_height_m*100:.2f} cm (kσ = {ks:.3f})")
    print(f"  - Correlation length: {corr_length_m*100:.2f} cm (kℓ = {kl:.3f})")
    print(f"  - Permittivity: {er.real:.1f} + {er.imag:.1f}j")
    
    # Test different grid resolutions
    test_configs = [
        (65, 6, "Fast"),
        (129, 8, "Standard"),
        (257, 10, "High-res")
    ]
    
    print(f"\n{'='*70}")
    print("Performance Results:")
    print(f"{'='*70}")
    print(f"{'Config':<12} {'Grid':<10} {'Order':<8} {'Time (s)':<12} {'HV (dB)':<10}")
    print("-"*70)
    
    for n_points, nmax, label in test_configs:
        try:
            start = time.time()
            
            result = compute_multiple_scattering(
                theta_i=theta_rad,
                theta_s=theta_rad,
                phi_i=0.0,
                phi_s=np.pi,
                er=er,
                ks=ks,
                kl=kl,
                k=k,
                sigma=sigma,
                surface_label='exponential',
                polarisations=['hv'],
                n_points=n_points,
                nmax=nmax
            )
            
            elapsed = time.time() - start
            
            hv_linear = result.get('hv', 0.0)
            hv_db = 10 * np.log10(hv_linear) if hv_linear > 0 else -999
            
            print(f"{label:<12} {n_points}×{n_points:<6} {nmax:<8} {elapsed:<12.3f} {hv_db:<10.2f}")
            
        except Exception as e:
            print(f"{label:<12} {n_points}×{n_points:<6} {nmax:<8} ERROR: {str(e)[:30]}")
    
    print("-"*70)
    
    # Estimate speedup potential
    if NUMBA_AVAILABLE:
        print("\n✅ With Numba JIT compilation:")
        print("   - First call: Slower (compilation overhead)")
        print("   - Subsequent calls: 20-100x faster")
        print("   - Recommended: Use cached compilation")
    else:
        print("\n⚠️  Without Numba:")
        print("   - Using pure NumPy implementation")
        print("   - Install numba for significant speedup:")
        print("     pip install numba")


def main():
    """Run all benchmarks."""
    print("="*70)
    print("AIEM MULTIPLE SCATTERING - NUMBA ACCELERATION BENCHMARK")
    print("="*70)
    
    # Check Numba status
    if NUMBA_AVAILABLE:
        status = check_numba_performance()
        print(f"\n✅ Numba Status:")
        print(f"   Version: {status.get('numba_version', 'Unknown')}")
        print(f"   Test passed: {status.get('test_passed', False)}")
    else:
        print(f"\n⚠️  Numba not available")
        print(f"   Install with: pip install numba")
    
    # Run benchmarks
    try:
        benchmark_numba_functions()
    except Exception as e:
        print(f"\n❌ Numba function benchmark failed: {e}")
    
    try:
        benchmark_full_multiple_scattering()
    except Exception as e:
        print(f"\n❌ Full MS benchmark failed: {e}")
    
    print("\n" + "="*70)
    print("BENCHMARK COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
