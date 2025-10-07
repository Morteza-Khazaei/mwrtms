"""
FINAL EMPIRICAL TEST: Compare AIEM and I2EM against NMM3D.
This will definitively show which model performs better for co-pol bands.
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
    """Load NMM3D reference data with CORRECT format."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 7:
                    theta = float(parts[0])
                    kL = float(parts[1])
                    eps_r = float(parts[2])
                    eps_i = float(parts[3])
                    ks = float(parts[4])
                    vv_db = float(parts[5])
                    hh_db = float(parts[6])
                    hv_db = float(parts[7]) if len(parts) > 7 and parts[7] != '-Inf' else None
                    
                    data.append({
                        'theta': theta,
                        'kL': kL,
                        'ks': ks,
                        'eps_r': eps_r,
                        'eps_i': eps_i,
                        'vv': vv_db,
                        'hh': hh_db,
                        'hv': hv_db
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
    corr = np.corrcoef(model_vals, nmm3d_vals)[0, 1] if len(model_vals) > 1 else np.nan
    
    return rmse, mae, bias, corr


def main():
    print("=" * 80)
    print("FINAL EMPIRICAL TEST: AIEM vs I2EM vs NMM3D")
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
    frequency_hz = 5.4e9  # C-band
    
    # Storage
    results_aiem = {'vv': [], 'hh': []}
    results_i2em = {'vv': [], 'hh': []}
    nmm3d_results = {'vv': [], 'hh': []}
    
    print("Running tests...")
    print()
    
    success_count = 0
    aiem_fail = 0
    i2em_fail = 0
    
    for idx, case in enumerate(nmm3d_data):
        theta_deg = case['theta']
        kL = case['kL']
        ks = case['ks']
        eps = complex(case['eps_r'], case['eps_i'])
        
        # Convert ks, kL to sigma, L
        k = 2 * np.pi * frequency_hz / 3e8
        sigma_m = ks / k
        L_m = kL / k
        
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
            
            # AIEM
            aiem_ok = False
            try:
                model_aiem = AIEMModel(
                    wave, geometry, surface,
                    correlation_type="exponential",
                    auto_terms=False,
                    spectral_terms=10,
                    use_i2em_transition=False
                )
                result_aiem = model_aiem.run(air, soil)
                aiem_vv = 10 * np.log10(result_aiem['vv'] + 1e-30)
                aiem_hh = 10 * np.log10(result_aiem['hh'] + 1e-30)
                if np.isfinite(aiem_vv) and np.isfinite(aiem_hh):
                    aiem_ok = True
            except:
                aiem_fail += 1
            
            # I2EM
            i2em_ok = False
            try:
                model_i2em = I2EMModel(
                    wave, geometry, surface,
                    correlation_type="exponential",
                    auto_terms=False,
                    spectral_terms=10
                )
                result_i2em = model_i2em.run(air, soil)
                i2em_vv = 10 * np.log10(result_i2em['vv'] + 1e-30)
                i2em_hh = 10 * np.log10(result_i2em['hh'] + 1e-30)
                if np.isfinite(i2em_vv) and np.isfinite(i2em_hh):
                    i2em_ok = True
            except:
                i2em_fail += 1
            
            # Store if both succeeded
            if aiem_ok and i2em_ok:
                results_aiem['vv'].append(aiem_vv)
                results_aiem['hh'].append(aiem_hh)
                results_i2em['vv'].append(i2em_vv)
                results_i2em['hh'].append(i2em_hh)
                nmm3d_results['vv'].append(case['vv'])
                nmm3d_results['hh'].append(case['hh'])
                success_count += 1
            
            if (idx + 1) % 50 == 0:
                print(f"  {idx+1}/{len(nmm3d_data)} processed, {success_count} successful")
        
        except Exception as e:
            continue
    
    print()
    print(f"✓ {success_count} cases where both models succeeded")
    print(f"  AIEM failures: {aiem_fail}")
    print(f"  I2EM failures: {i2em_fail}")
    print()
    
    if success_count < 10:
        print("✗ Not enough data")
        return
    
    # VV Results
    print("=" * 80)
    print("VV POLARIZATION RESULTS")
    print("=" * 80)
    print()
    
    rmse_aiem_vv, mae_aiem_vv, bias_aiem_vv, corr_aiem_vv = compute_metrics(
        results_aiem['vv'], nmm3d_results['vv']
    )
    rmse_i2em_vv, mae_i2em_vv, bias_i2em_vv, corr_i2em_vv = compute_metrics(
        results_i2em['vv'], nmm3d_results['vv']
    )
    
    print(f"AIEM:  RMSE={rmse_aiem_vv:6.2f} dB  Bias={bias_aiem_vv:+6.2f} dB  Corr={corr_aiem_vv:.3f}")
    print(f"I2EM:  RMSE={rmse_i2em_vv:6.2f} dB  Bias={bias_i2em_vv:+6.2f} dB  Corr={corr_i2em_vv:.3f}")
    print()
    print(f"Δ(AIEM-I2EM): RMSE={rmse_aiem_vv-rmse_i2em_vv:+6.2f} dB  Bias={bias_aiem_vv-bias_i2em_vv:+6.2f} dB")
    print()
    
    if rmse_i2em_vv < rmse_aiem_vv:
        print(f"✓ I2EM WINS for VV (RMSE {rmse_aiem_vv-rmse_i2em_vv:.2f} dB lower)")
    else:
        print(f"✗ AIEM wins for VV")
    print()
    
    # HH Results
    print("=" * 80)
    print("HH POLARIZATION RESULTS")
    print("=" * 80)
    print()
    
    rmse_aiem_hh, mae_aiem_hh, bias_aiem_hh, corr_aiem_hh = compute_metrics(
        results_aiem['hh'], nmm3d_results['hh']
    )
    rmse_i2em_hh, mae_i2em_hh, bias_i2em_hh, corr_i2em_hh = compute_metrics(
        results_i2em['hh'], nmm3d_results['hh']
    )
    
    print(f"AIEM:  RMSE={rmse_aiem_hh:6.2f} dB  Bias={bias_aiem_hh:+6.2f} dB  Corr={corr_aiem_hh:.3f}")
    print(f"I2EM:  RMSE={rmse_i2em_hh:6.2f} dB  Bias={bias_i2em_hh:+6.2f} dB  Corr={corr_i2em_hh:.3f}")
    print()
    print(f"Δ(AIEM-I2EM): RMSE={rmse_aiem_hh-rmse_i2em_hh:+6.2f} dB  Bias={bias_aiem_hh-bias_i2em_hh:+6.2f} dB")
    print()
    
    if rmse_i2em_hh < rmse_aiem_hh:
        print(f"✓ I2EM WINS for HH (RMSE {rmse_aiem_hh-rmse_i2em_hh:.2f} dB lower)")
    else:
        print(f"✗ AIEM wins for HH")
    print()
    
    # Final verdict
    print("=" * 80)
    print("EMPIRICAL VERDICT")
    print("=" * 80)
    print()
    
    if rmse_i2em_vv < rmse_aiem_vv and rmse_i2em_hh < rmse_aiem_hh:
        print("✓✓ CONFIRMED: I2EM performs better for BOTH co-pol bands!")
        print()
        print("Why?")
        if bias_aiem_vv > 1.0:
            print(f"  • AIEM has positive bias: VV={bias_aiem_vv:+.2f} dB, HH={bias_aiem_hh:+.2f} dB")
            print("  • This confirms AIEM's transition function causes over-prediction")
        print(f"  • I2EM has lower bias: VV={bias_i2em_vv:+.2f} dB, HH={bias_i2em_hh:+.2f} dB")
        print()
        print("Solution: Use I2EM transition in AIEM")
        print("  model = AIEMModel(..., use_i2em_transition=True)")
    else:
        print("Results show mixed or unexpected performance")
    
    print()


if __name__ == "__main__":
    main()
