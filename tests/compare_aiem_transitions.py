"""
Compare AIEM with original transition vs I2EM-style transition.

This script tests whether using the I2EM transition method in AIEM
reduces the systematic bias in co-polarization channels.
"""

import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mwrtms.core import ElectromagneticWave, ScatteringGeometry, PolarizationState
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


def compute_metrics(aiem_vals, nmm3d_vals):
    """Compute RMSE, MAE, and Bias."""
    aiem_vals = np.array(aiem_vals)
    nmm3d_vals = np.array(nmm3d_vals)
    
    diff = aiem_vals - nmm3d_vals
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))
    bias = np.mean(diff)
    corr = np.corrcoef(aiem_vals, nmm3d_vals)[0, 1]
    
    return rmse, mae, bias, corr


def main():
    print("=" * 80)
    print("AIEM Transition Method Comparison")
    print("=" * 80)
    print()
    
    # Load NMM3D data
    nmm3d_data = load_nmm3d_data()
    print(f"Loaded {len(nmm3d_data)} NMM3D reference points")
    print()
    
    # Test parameters
    frequency_hz = 5.4e9  # 5.4 GHz (C-band)
    theta_deg = 40.0
    
    # Storage for results
    results_original = {'vv': [], 'hh': [], 'hv': []}
    results_i2em_trans = {'vv': [], 'hh': [], 'hv': []}
    results_i2em_model = {'vv': [], 'hh': [], 'hv': []}
    nmm3d_results = {'vv': [], 'hh': [], 'hv': []}
    
    print("Computing backscatter for all test cases...")
    print()
    
    for idx, case in enumerate(nmm3d_data[:20]):  # Test first 20 cases
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
            
            # Test 1: AIEM with original transition
            model_orig = AIEMModel(
                wave, geometry, surface,
                correlation_type="exponential",
                auto_terms=False,
                spectral_terms=15,
                use_i2em_transition=False
            )
            result_orig = model_orig.run(air, soil)
            
            # Test 2: AIEM with I2EM transition
            model_i2em_trans = AIEMModel(
                wave, geometry, surface,
                correlation_type="exponential",
                auto_terms=False,
                spectral_terms=15,
                use_i2em_transition=True
            )
            result_i2em_trans = model_i2em_trans.run(air, soil)
            
            # Test 3: I2EM model for reference
            model_i2em = I2EMModel(
                wave, geometry, surface,
                correlation_type="exponential",
                auto_terms=False,
                spectral_terms=15
            )
            result_i2em = model_i2em.run(air, soil)
            
            # Store results (convert to dB)
            result_orig_db = result_orig.to_db()
            results_original['vv'].append(result_orig_db['vv'])
            results_original['hh'].append(result_orig_db['hh'])
            results_original['hv'].append(result_orig_db['hv'])
            
            result_i2em_trans_db = result_i2em_trans.to_db()
            results_i2em_trans['vv'].append(result_i2em_trans_db['vv'])
            results_i2em_trans['hh'].append(result_i2em_trans_db['hh'])
            results_i2em_trans['hv'].append(result_i2em_trans_db['hv'])
            
            result_i2em_db = result_i2em.to_db()
            results_i2em_model['vv'].append(result_i2em_db['vv'])
            results_i2em_model['hh'].append(result_i2em_db['hh'])
            results_i2em_model['hv'].append(result_i2em_db['hv'])
            
            nmm3d_results['vv'].append(case['vv'])
            nmm3d_results['hh'].append(case['hh'])
            if case['hv'] is not None and not np.isnan(case['hv']) and np.isfinite(case['hv']):
                nmm3d_results['hv'].append(case['hv'])
            
            if (idx + 1) % 5 == 0:
                print(f"  Processed {idx + 1}/{min(20, len(nmm3d_data))} cases...")
        except Exception as e:
            print(f"  Skipping case {idx + 1} due to error: {e}")
            continue
    
    print()
    print("=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print()
    
    # VV Polarization
    print("VV POLARIZATION:")
    print("-" * 80)
    
    rmse_orig, mae_orig, bias_orig, corr_orig = compute_metrics(
        results_original['vv'], nmm3d_results['vv']
    )
    print(f"AIEM (Original Transition):")
    print(f"  RMSE: {rmse_orig:.2f} dB  |  MAE: {mae_orig:.2f} dB  |  Bias: {bias_orig:+.2f} dB  |  Corr: {corr_orig:.3f}")
    
    rmse_i2em, mae_i2em, bias_i2em, corr_i2em = compute_metrics(
        results_i2em_trans['vv'], nmm3d_results['vv']
    )
    print(f"AIEM (I2EM Transition):     ")
    print(f"  RMSE: {rmse_i2em:.2f} dB  |  MAE: {mae_i2em:.2f} dB  |  Bias: {bias_i2em:+.2f} dB  |  Corr: {corr_i2em:.3f}")
    
    rmse_i2em_model, mae_i2em_model, bias_i2em_model, corr_i2em_model = compute_metrics(
        results_i2em_model['vv'], nmm3d_results['vv']
    )
    print(f"I2EM Model (Reference):     ")
    print(f"  RMSE: {rmse_i2em_model:.2f} dB  |  MAE: {mae_i2em_model:.2f} dB  |  Bias: {bias_i2em_model:+.2f} dB  |  Corr: {corr_i2em_model:.3f}")
    
    print()
    print(f"VV Improvement (RMSE): {rmse_orig - rmse_i2em:+.2f} dB")
    print(f"VV Improvement (Bias): {abs(bias_orig) - abs(bias_i2em):+.2f} dB")
    print()
    
    # HH Polarization
    print("HH POLARIZATION:")
    print("-" * 80)
    
    rmse_orig, mae_orig, bias_orig, corr_orig = compute_metrics(
        results_original['hh'], nmm3d_results['hh']
    )
    print(f"AIEM (Original Transition):")
    print(f"  RMSE: {rmse_orig:.2f} dB  |  MAE: {mae_orig:.2f} dB  |  Bias: {bias_orig:+.2f} dB  |  Corr: {corr_orig:.3f}")
    
    rmse_i2em, mae_i2em, bias_i2em, corr_i2em = compute_metrics(
        results_i2em_trans['hh'], nmm3d_results['hh']
    )
    print(f"AIEM (I2EM Transition):     ")
    print(f"  RMSE: {rmse_i2em:.2f} dB  |  MAE: {mae_i2em:.2f} dB  |  Bias: {bias_i2em:+.2f} dB  |  Corr: {corr_i2em:.3f}")
    
    rmse_i2em_model, mae_i2em_model, bias_i2em_model, corr_i2em_model = compute_metrics(
        results_i2em_model['hh'], nmm3d_results['hh']
    )
    print(f"I2EM Model (Reference):     ")
    print(f"  RMSE: {rmse_i2em_model:.2f} dB  |  MAE: {mae_i2em_model:.2f} dB  |  Bias: {bias_i2em_model:+.2f} dB  |  Corr: {corr_i2em_model:.3f}")
    
    print()
    print(f"HH Improvement (RMSE): {rmse_orig - rmse_i2em:+.2f} dB")
    print(f"HH Improvement (Bias): {abs(bias_orig) - abs(bias_i2em):+.2f} dB")
    print()
    
    print("=" * 80)
    print("CONCLUSION:")
    print("=" * 80)
    if rmse_i2em < rmse_orig:
        print("✓ I2EM transition method IMPROVES AIEM co-pol accuracy!")
    else:
        print("✗ I2EM transition method does NOT improve AIEM co-pol accuracy.")
    print()


if __name__ == "__main__":
    main()
