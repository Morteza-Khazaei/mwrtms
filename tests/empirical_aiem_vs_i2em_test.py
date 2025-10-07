"""
Empirical test: Compare AIEM and I2EM against NMM3D reference data.
This test will show WHY I2EM performs better for co-pol bands.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mwrtms.core import ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface.iem.aiem import AIEMModel
from mwrtms.scattering.surface.iem.i2em import I2EMModel


def load_nmm3d_data(filepath="data/NMM3D_LUT_NRCS_40degree.dat"):
    """Load NMM3D reference data."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 6:
                    data.append({
                        'sigma': float(parts[0]),
                        'L': float(parts[1]),
                        'eps_r': float(parts[2]),
                        'eps_i': float(parts[3]),
                        'vv': float(parts[4]),
                        'hh': float(parts[5]),
                        'hv': float(parts[6]) if len(parts) > 6 else None
                    })
    return data


def compute_metrics(model_vals, nmm3d_vals):
    """Compute RMSE, MAE, Bias, and Correlation."""
    model_vals = np.array(model_vals)
    nmm3d_vals = np.array(nmm3d_vals)
    
    diff = model_vals - nmm3d_vals
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))
    bias = np.mean(diff)
    
    if len(model_vals) > 1:
        corr = np.corrcoef(model_vals, nmm3d_vals)[0, 1]
    else:
        corr = np.nan
    
    return rmse, mae, bias, corr


def main():
    print("=" * 80)
    print("EMPIRICAL TEST: AIEM vs I2EM vs NMM3D")
    print("=" * 80)
    print()
    print("This test will run both models against NMM3D reference data")
    print("to empirically demonstrate why I2EM performs better for co-pol bands.")
    print()
    
    # Load NMM3D data
    try:
        nmm3d_data = load_nmm3d_data()
        print(f"✓ Loaded {len(nmm3d_data)} NMM3D reference points")
    except FileNotFoundError:
        print("✗ NMM3D data file not found!")
        print("  Expected: data/NMM3D_LUT_NRCS_40degree.dat")
        return
    
    print()
    print("Test Configuration:")
    print("  Frequency: 5.4 GHz (C-band)")
    print("  Incidence angle: 40°")
    print("  Spectral terms: 15 (fixed)")
    print("  Correlation: exponential")
    print()
    
    # Test parameters
    frequency_hz = 5.4e9
    theta_deg = 40.0
    
    # Storage for results
    results_aiem = {'vv': [], 'hh': []}
    results_i2em = {'vv': [], 'hh': []}
    nmm3d_results = {'vv': [], 'hh': []}
    
    # Test parameters for detailed analysis
    test_cases = []
    
    print("Running tests on subset of NMM3D data...")
    print()
    
    # Test on first 30 cases for speed
    n_test = min(30, len(nmm3d_data))
    success_count = 0
    
    for idx, case in enumerate(nmm3d_data[:n_test]):
        sigma_m = case['sigma']
        L_m = case['L']
        eps = complex(case['eps_r'], case['eps_i'])
        
        # Skip invalid cases
        if sigma_m <= 0 or L_m <= 0:
            continue
        
        try:
            # Create wave and geometry
            wave = ElectromagneticWave(frequency_hz)
            geometry = ScatteringGeometry(theta_i_deg=theta_deg)
            
            # Create surface
            surface = build_surface_from_statistics(
                rms_height_m=sigma_m,
                correlation_length_m=L_m,
                correlation_type="exponential"
            )
            
            # Create media
            air = HomogeneousMedium(1.0 + 0.0j)
            soil = HomogeneousMedium(eps)
            
            # Test AIEM (original transition)
            model_aiem = AIEMModel(
                wave, geometry, surface,
                correlation_type="exponential",
                auto_terms=False,
                spectral_terms=15,
                use_i2em_transition=False  # Original AIEM transition
            )
            result_aiem = model_aiem.run(air, soil)
            result_aiem_db = result_aiem.to_db()
            
            # Test I2EM
            model_i2em = I2EMModel(
                wave, geometry, surface,
                correlation_type="exponential",
                auto_terms=False,
                spectral_terms=15
            )
            result_i2em = model_i2em.run(air, soil)
            result_i2em_db = result_i2em.to_db()
            
            # Store results
            results_aiem['vv'].append(result_aiem_db['vv'])
            results_aiem['hh'].append(result_aiem_db['hh'])
            
            results_i2em['vv'].append(result_i2em_db['vv'])
            results_i2em['hh'].append(result_i2em_db['hh'])
            
            nmm3d_results['vv'].append(case['vv'])
            nmm3d_results['hh'].append(case['hh'])
            
            # Store case details for analysis
            test_cases.append({
                'idx': idx,
                'sigma': sigma_m,
                'L': L_m,
                'eps': eps,
                'ks': 2 * np.pi * frequency_hz / 3e8 * sigma_m,
                'kL': 2 * np.pi * frequency_hz / 3e8 * L_m,
                'aiem_vv': result_aiem_db['vv'],
                'aiem_hh': result_aiem_db['hh'],
                'i2em_vv': result_i2em_db['vv'],
                'i2em_hh': result_i2em_db['hh'],
                'nmm3d_vv': case['vv'],
                'nmm3d_hh': case['hh'],
            })
            
            success_count += 1
            
            if (idx + 1) % 10 == 0:
                print(f"  Processed {idx + 1}/{n_test} cases... ({success_count} successful)")
        
        except Exception as e:
            print(f"  ✗ Case {idx + 1} failed: {str(e)[:50]}")
            continue
    
    print()
    print(f"✓ Successfully processed {success_count} cases")
    print()
    
    if success_count < 5:
        print("✗ Not enough successful cases for meaningful analysis")
        return
    
    # Compute metrics
    print("=" * 80)
    print("RESULTS: VV POLARIZATION")
    print("=" * 80)
    print()
    
    rmse_aiem_vv, mae_aiem_vv, bias_aiem_vv, corr_aiem_vv = compute_metrics(
        results_aiem['vv'], nmm3d_results['vv']
    )
    rmse_i2em_vv, mae_i2em_vv, bias_i2em_vv, corr_i2em_vv = compute_metrics(
        results_i2em['vv'], nmm3d_results['vv']
    )
    
    print(f"AIEM (Original Transition):")
    print(f"  RMSE: {rmse_aiem_vv:7.2f} dB")
    print(f"  MAE:  {mae_aiem_vv:7.2f} dB")
    print(f"  Bias: {bias_aiem_vv:+7.2f} dB")
    print(f"  Corr: {corr_aiem_vv:7.3f}")
    print()
    
    print(f"I2EM:")
    print(f"  RMSE: {rmse_i2em_vv:7.2f} dB")
    print(f"  MAE:  {mae_i2em_vv:7.2f} dB")
    print(f"  Bias: {bias_i2em_vv:+7.2f} dB")
    print(f"  Corr: {corr_i2em_vv:7.3f}")
    print()
    
    print(f"Difference (AIEM - I2EM):")
    print(f"  ΔRMSE: {rmse_aiem_vv - rmse_i2em_vv:+7.2f} dB")
    print(f"  ΔBias: {bias_aiem_vv - bias_i2em_vv:+7.2f} dB")
    print()
    
    if abs(bias_i2em_vv) < abs(bias_aiem_vv):
        print(f"✓ I2EM has LOWER bias ({abs(bias_i2em_vv):.2f} vs {abs(bias_aiem_vv):.2f} dB)")
    else:
        print(f"✗ AIEM has lower bias")
    print()
    
    # HH Polarization
    print("=" * 80)
    print("RESULTS: HH POLARIZATION")
    print("=" * 80)
    print()
    
    rmse_aiem_hh, mae_aiem_hh, bias_aiem_hh, corr_aiem_hh = compute_metrics(
        results_aiem['hh'], nmm3d_results['hh']
    )
    rmse_i2em_hh, mae_i2em_hh, bias_i2em_hh, corr_i2em_hh = compute_metrics(
        results_i2em['hh'], nmm3d_results['hh']
    )
    
    print(f"AIEM (Original Transition):")
    print(f"  RMSE: {rmse_aiem_hh:7.2f} dB")
    print(f"  MAE:  {mae_aiem_hh:7.2f} dB")
    print(f"  Bias: {bias_aiem_hh:+7.2f} dB")
    print(f"  Corr: {corr_aiem_hh:7.3f}")
    print()
    
    print(f"I2EM:")
    print(f"  RMSE: {rmse_i2em_hh:7.2f} dB")
    print(f"  MAE:  {mae_i2em_hh:7.2f} dB")
    print(f"  Bias: {bias_i2em_hh:+7.2f} dB")
    print(f"  Corr: {corr_i2em_hh:7.3f}")
    print()
    
    print(f"Difference (AIEM - I2EM):")
    print(f"  ΔRMSE: {rmse_aiem_hh - rmse_i2em_hh:+7.2f} dB")
    print(f"  ΔBias: {bias_aiem_hh - bias_i2em_hh:+7.2f} dB")
    print()
    
    if abs(bias_i2em_hh) < abs(bias_aiem_hh):
        print(f"✓ I2EM has LOWER bias ({abs(bias_i2em_hh):.2f} vs {abs(bias_aiem_hh):.2f} dB)")
    else:
        print(f"✗ AIEM has lower bias")
    print()
    
    # Detailed case analysis
    print("=" * 80)
    print("DETAILED CASE ANALYSIS (First 5 cases)")
    print("=" * 80)
    print()
    
    for i, tc in enumerate(test_cases[:5]):
        print(f"Case {tc['idx'] + 1}: σ={tc['sigma']*100:.1f}cm, L={tc['L']*100:.1f}cm, "
              f"ks={tc['ks']:.2f}, kL={tc['kL']:.1f}")
        print(f"  VV: NMM3D={tc['nmm3d_vv']:7.2f} dB  |  "
              f"AIEM={tc['aiem_vv']:7.2f} dB (Δ={tc['aiem_vv']-tc['nmm3d_vv']:+6.2f})  |  "
              f"I2EM={tc['i2em_vv']:7.2f} dB (Δ={tc['i2em_vv']-tc['nmm3d_vv']:+6.2f})")
        print(f"  HH: NMM3D={tc['nmm3d_hh']:7.2f} dB  |  "
              f"AIEM={tc['aiem_hh']:7.2f} dB (Δ={tc['aiem_hh']-tc['nmm3d_hh']:+6.2f})  |  "
              f"I2EM={tc['i2em_hh']:7.2f} dB (Δ={tc['i2em_hh']-tc['nmm3d_hh']:+6.2f})")
        print()
    
    # Summary
    print("=" * 80)
    print("EMPIRICAL CONCLUSION")
    print("=" * 80)
    print()
    
    vv_winner = "I2EM" if abs(bias_i2em_vv) < abs(bias_aiem_vv) else "AIEM"
    hh_winner = "I2EM" if abs(bias_i2em_hh) < abs(bias_aiem_hh) else "AIEM"
    
    print(f"VV Polarization: {vv_winner} performs better")
    print(f"  - Lower bias: {min(abs(bias_i2em_vv), abs(bias_aiem_vv)):.2f} dB vs "
          f"{max(abs(bias_i2em_vv), abs(bias_aiem_vv)):.2f} dB")
    print()
    
    print(f"HH Polarization: {hh_winner} performs better")
    print(f"  - Lower bias: {min(abs(bias_i2em_hh), abs(bias_aiem_hh)):.2f} dB vs "
          f"{max(abs(bias_i2em_hh), abs(bias_aiem_hh)):.2f} dB")
    print()
    
    if vv_winner == "I2EM" and hh_winner == "I2EM":
        print("✓ EMPIRICALLY CONFIRMED: I2EM performs better for BOTH co-pol bands")
        print()
        print("Why? Based on the data:")
        if bias_aiem_vv > 0 and bias_aiem_hh > 0:
            print("  • AIEM shows POSITIVE bias (over-prediction)")
            print("  • I2EM shows bias closer to zero")
            print("  • This confirms AIEM's transition function is too aggressive")
        print()
        print("The solution: Use I2EM transition in AIEM")
        print("  model = AIEMModel(..., use_i2em_transition=True)")
    else:
        print("⚠ Results differ from expected - further investigation needed")
    
    print()
    print("=" * 80)


if __name__ == "__main__":
    main()
