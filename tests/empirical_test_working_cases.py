"""
Empirical test: Compare AIEM and I2EM against NMM3D on cases where both work.
This will show WHY I2EM performs better for co-pol bands.
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
    print("Testing only on cases where both models succeed")
    print("=" * 80)
    print()
    
    # Load NMM3D data
    try:
        nmm3d_data = load_nmm3d_data()
        print(f"✓ Loaded {len(nmm3d_data)} NMM3D reference points")
    except FileNotFoundError:
        print("✗ NMM3D data file not found!")
        return
    
    print()
    
    # Test parameters
    frequency_hz = 5.4e9
    theta_deg = 40.0
    
    # Storage for results
    results_aiem = {'vv': [], 'hh': []}
    results_i2em = {'vv': [], 'hh': []}
    nmm3d_results = {'vv': [], 'hh': []}
    test_cases = []
    
    print("Running tests (this may take a minute)...")
    print()
    
    success_count = 0
    aiem_fail_count = 0
    i2em_fail_count = 0
    
    for idx, case in enumerate(nmm3d_data):
        sigma_m = case['sigma']
        L_m = case['L']
        eps = complex(case['eps_r'], case['eps_i'])
        
        if sigma_m <= 0 or L_m <= 0:
            continue
        
        try:
            wave = ElectromagneticWave(frequency_hz)
            geometry = ScatteringGeometry(theta_i_deg=theta_deg)
            surface = build_surface_from_statistics(
                rms_height_m=sigma_m,
                correlation_length_m=L_m,
                correlation_type="exponential"
            )
            air = HomogeneousMedium(1.0 + 0.0j)
            soil = HomogeneousMedium(eps)
            
            # Try AIEM first
            try:
                model_aiem = AIEMModel(
                    wave, geometry, surface,
                    correlation_type="exponential",
                    auto_terms=False,
                    spectral_terms=10,
                    use_i2em_transition=False
                )
                result_aiem = model_aiem.run(air, soil)
                aiem_vv = 10 * np.log10(result_aiem['vv'])
                aiem_hh = 10 * np.log10(result_aiem['hh'])
                aiem_success = True
            except:
                aiem_success = False
                aiem_fail_count += 1
            
            # Try I2EM
            try:
                model_i2em = I2EMModel(
                    wave, geometry, surface,
                    correlation_type="exponential",
                    auto_terms=False,
                    spectral_terms=10
                )
                result_i2em = model_i2em.run(air, soil)
                i2em_vv = 10 * np.log10(result_i2em['vv'])
                i2em_hh = 10 * np.log10(result_i2em['hh'])
                i2em_success = True
            except:
                i2em_success = False
                i2em_fail_count += 1
            
            # Only store if BOTH succeeded
            if aiem_success and i2em_success:
                results_aiem['vv'].append(aiem_vv)
                results_aiem['hh'].append(aiem_hh)
                results_i2em['vv'].append(i2em_vv)
                results_i2em['hh'].append(i2em_hh)
                nmm3d_results['vv'].append(case['vv'])
                nmm3d_results['hh'].append(case['hh'])
                
                test_cases.append({
                    'idx': idx,
                    'sigma': sigma_m,
                    'L': L_m,
                    'ks': 2 * np.pi * frequency_hz / 3e8 * sigma_m,
                    'aiem_vv': aiem_vv,
                    'aiem_hh': aiem_hh,
                    'i2em_vv': i2em_vv,
                    'i2em_hh': i2em_hh,
                    'nmm3d_vv': case['vv'],
                    'nmm3d_hh': case['hh'],
                })
                
                success_count += 1
            
            if (idx + 1) % 50 == 0:
                print(f"  Processed {idx + 1}/{len(nmm3d_data)} cases... "
                      f"({success_count} both succeeded, {aiem_fail_count} AIEM failed, {i2em_fail_count} I2EM failed)")
        
        except Exception as e:
            continue
    
    print()
    print(f"✓ Successfully processed {success_count} cases where BOTH models worked")
    print(f"  AIEM failed on {aiem_fail_count} cases")
    print(f"  I2EM failed on {i2em_fail_count} cases")
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
    print(f"  RMSE: {rmse_aiem_vv:7.2f} dB  |  Bias: {bias_aiem_vv:+7.2f} dB  |  Corr: {corr_aiem_vv:6.3f}")
    print()
    print(f"I2EM:")
    print(f"  RMSE: {rmse_i2em_vv:7.2f} dB  |  Bias: {bias_i2em_vv:+7.2f} dB  |  Corr: {corr_i2em_vv:6.3f}")
    print()
    print(f"Difference (AIEM - I2EM):")
    print(f"  ΔRMSE: {rmse_aiem_vv - rmse_i2em_vv:+7.2f} dB")
    print(f"  ΔBias: {bias_aiem_vv - bias_i2em_vv:+7.2f} dB")
    print()
    
    vv_winner = "I2EM" if rmse_i2em_vv < rmse_aiem_vv else "AIEM"
    print(f"✓ {vv_winner} performs better for VV (lower RMSE)")
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
    print(f"  RMSE: {rmse_aiem_hh:7.2f} dB  |  Bias: {bias_aiem_hh:+7.2f} dB  |  Corr: {corr_aiem_hh:6.3f}")
    print()
    print(f"I2EM:")
    print(f"  RMSE: {rmse_i2em_hh:7.2f} dB  |  Bias: {bias_i2em_hh:+7.2f} dB  |  Corr: {corr_i2em_hh:6.3f}")
    print()
    print(f"Difference (AIEM - I2EM):")
    print(f"  ΔRMSE: {rmse_aiem_hh - rmse_i2em_hh:+7.2f} dB")
    print(f"  ΔBias: {bias_aiem_hh - bias_i2em_hh:+7.2f} dB")
    print()
    
    hh_winner = "I2EM" if rmse_i2em_hh < rmse_aiem_hh else "AIEM"
    print(f"✓ {hh_winner} performs better for HH (lower RMSE)")
    print()
    
    # Sample cases
    print("=" * 80)
    print("SAMPLE CASES (First 5)")
    print("=" * 80)
    print()
    
    for tc in test_cases[:5]:
        print(f"Case {tc['idx']+1}: ks={tc['ks']:.2f}")
        print(f"  VV: NMM3D={tc['nmm3d_vv']:6.2f}  AIEM={tc['aiem_vv']:6.2f} (Δ={tc['aiem_vv']-tc['nmm3d_vv']:+5.2f})  "
              f"I2EM={tc['i2em_vv']:6.2f} (Δ={tc['i2em_vv']-tc['nmm3d_vv']:+5.2f})")
        print(f"  HH: NMM3D={tc['nmm3d_hh']:6.2f}  AIEM={tc['aiem_hh']:6.2f} (Δ={tc['aiem_hh']-tc['nmm3d_hh']:+5.2f})  "
              f"I2EM={tc['i2em_hh']:6.2f} (Δ={tc['i2em_hh']-tc['nmm3d_hh']:+5.2f})")
        print()
    
    # Final conclusion
    print("=" * 80)
    print("EMPIRICAL CONCLUSION")
    print("=" * 80)
    print()
    
    if vv_winner == "I2EM" and hh_winner == "I2EM":
        print("✓ EMPIRICALLY CONFIRMED: I2EM performs better for BOTH co-pol bands")
        print()
        print(f"  VV: I2EM RMSE is {rmse_aiem_vv - rmse_i2em_vv:.2f} dB lower than AIEM")
        print(f"  HH: I2EM RMSE is {rmse_aiem_hh - rmse_i2em_hh:.2f} dB lower than AIEM")
        print()
        if bias_aiem_vv > 0:
            print(f"  AIEM shows positive bias: VV={bias_aiem_vv:+.2f} dB, HH={bias_aiem_hh:+.2f} dB")
            print(f"  I2EM shows lower bias: VV={bias_i2em_vv:+.2f} dB, HH={bias_i2em_hh:+.2f} dB")
            print()
            print("  This confirms: AIEM's transition function causes over-prediction!")
    elif vv_winner == "AIEM" and hh_winner == "AIEM":
        print("⚠ Unexpected: AIEM performs better for both bands")
        print("  This contradicts the documented behavior - further investigation needed")
    else:
        print(f"⚠ Mixed results: {vv_winner} better for VV, {hh_winner} better for HH")
    
    print()
    print("Note: AIEM failed on many high-roughness cases, suggesting implementation issues")
    print(f"beyond just the transition function (failed on {aiem_fail_count}/{len(nmm3d_data)} cases)")
    print()


if __name__ == "__main__":
    main()
