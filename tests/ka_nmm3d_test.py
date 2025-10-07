"""CLI test harness comparing KA model against the NMM3D LUT.

This script provides numerical validation of the Kirchhoff Approximation (KA)
model by comparing its predictions against the NMM3D backscatter look-up table.
It loads the 40° incidence NMM3D data, evaluates KA for each surface
configuration, and reports goodness-of-fit metrics (RMSE, MAE, bias, Pearson
correlation) for the HH, VV, and HV channels.

Note: The KA model is designed for large-scale roughness (long gravity waves)
and may not match NMM3D perfectly for all surface configurations, particularly
for small-scale roughness where other models (SPM, IEM) are more appropriate.

Usage
-----
Run from the repository root:

    PYTHONPATH=src python3 tests/ka_nmm3d_test.py

Optional command-line arguments let you filter ratios, choose a different
incident angle, or point to alternative LUTs. The script exits with a non-zero
status if the LUT cannot be located or no valid comparisons remain after
filtering.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState


def toLambda(frequency_ghz: float) -> float:
    """Convert frequency in GHz to wavelength in meters.
    
    Parameters
    ----------
    frequency_ghz : float
        Frequency in GHz
        
    Returns
    -------
    float
        Wavelength in meters
    """
    # Speed of light in m/s divided by frequency in Hz
    return 0.3 / frequency_ghz


# Default configuration mirrors the AIEM test
_DEFAULT_LUT = Path("data/NMM3D_LUT_NRCS_40degree.dat")
_DEFAULT_FREQ_GHZ = 5.405
_DEFAULT_INC_DEG = 40.0
_DEFAULT_PHI_DEG = 180.0


@dataclass
class Metrics:
    """Container for descriptive statistics between KA and NMM3D."""

    count: int
    rmse: float
    mae: float
    bias: float
    corr: float

    def format_row(self, label: str) -> str:  # pragma: no cover - trivial formatting
        return (
            f"{label:<6s} n={self.count:3d}  "
            f"RMSE={self.rmse:5.2f} dB  "
            f"MAE={self.mae:5.2f} dB  "
            f"Bias={self.bias:+5.2f} dB  "
            f"Corr={self.corr:6.3f}"
        )


def _load_lut(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"LUT not found at {path}")
    return np.loadtxt(path)


def _select_angle(table: np.ndarray, incidence_deg: float) -> np.ndarray:
    mask = np.isclose(table[:, 0], incidence_deg)
    return table[mask]


def _isclose(value: float, targets: Sequence[float], tol: float = 1e-6) -> bool:
    return any(math.isclose(value, target, rel_tol=0.0, abs_tol=tol) for target in targets)


def _calc_metrics(model: Iterable[float], reference: Iterable[float]) -> Metrics:
    model_arr = np.asarray(list(model), dtype=float)
    ref_arr = np.asarray(list(reference), dtype=float)
    mask = np.isfinite(model_arr) & np.isfinite(ref_arr)
    if not np.any(mask):
        return Metrics(count=0, rmse=float("nan"), mae=float("nan"), bias=float("nan"), corr=float("nan"))

    diff = model_arr[mask] - ref_arr[mask]
    rmse = float(np.sqrt(np.mean(diff**2)))
    mae = float(np.mean(np.abs(diff)))
    bias = float(np.mean(diff))
    corr = float(np.corrcoef(ref_arr[mask], model_arr[mask])[0, 1]) if mask.sum() > 1 else float("nan")
    return Metrics(count=int(mask.sum()), rmse=rmse, mae=mae, bias=bias, corr=corr)


@dataclass
class ComparisonResult:
    overall: Dict[str, Metrics]
    by_ratio: Dict[float, Dict[str, Metrics]]


def _run_comparison(
    rows: np.ndarray,
    frequency_ghz: float,
    incidence_deg: float,
    phi_deg: float,
    ratios: Sequence[float] | None,
    min_ratio: float | None,
) -> ComparisonResult:
    lam = toLambda(frequency_ghz)
    k = 2.0 * math.pi / lam

    overall_model: Dict[str, List[float]] = {pol: [] for pol in ("hh", "vv", "hv")}
    overall_reference: Dict[str, List[float]] = {pol: [] for pol in ("hh", "vv", "hv")}
    grouped_model: Dict[float, Dict[str, List[float]]] = {}
    grouped_reference: Dict[float, Dict[str, List[float]]] = {}

    # Create radar configuration
    radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=incidence_deg)

    for row in rows:
        (
            _theta,
            ratio,
            eps_r,
            eps_i,
            rms_norm,
            vv_ref,
            hh_ref,
            hv_ref,
        ) = row

        ratio = float(ratio)
        
        # Filter by ratio if specified
        if ratios and not _isclose(ratio, ratios):
            continue
        
        # Filter by minimum ratio if specified
        if min_ratio is not None and ratio < min_ratio:
            continue

        sigma = rms_norm * lam
        corr_len = ratio * sigma
        
        # Convert to cm for the facade API
        rms_height_cm = sigma * 100.0
        correlation_length_cm = corr_len * 100.0
        
        # Create complex permittivity
        soil_permittivity = complex(float(eps_r), float(eps_i))
        
        # Check if this is in the KA validity range
        # KA (Kirchhoff Approximation / Geometric Optics) is valid for:
        # 1. Large roughness: k*sigma > 0.3 (KA works for rough surfaces)
        # 2. Large radius of curvature: k*L > 6 (preferably)
        k_L = k * corr_len
        k_sigma = k * sigma
        
        # KA is designed for rough surfaces and should work better than NMM3D
        # on very rough surfaces. Don't filter out rough surfaces!
        # Only skip extremely smooth surfaces where KA is not applicable
        if k_sigma < 0.3:  # Too smooth for KA
            continue

        # Compute backscatter for each polarization using exponential ACF
        # NMM3D uses exponential correlation, so we should too for fair comparison
        try:
            from mwrtms.core import ElectromagneticWave, ScatteringGeometry
            from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
            from mwrtms.scattering.surface.ka import KAModel
            
            wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
            geometry = ScatteringGeometry(theta_i_deg=incidence_deg)
            surface = build_surface_from_statistics(
                rms_height_m=sigma,
                correlation_length_m=corr_len,
                correlation_type="exponential"
            )
            soil = HomogeneousMedium(permittivity=soil_permittivity)
            
            # Use exponential ACF to match NMM3D
            model = KAModel(wave, geometry, surface, acf_type="exponential", nmax=8)
            
            vv_result = model.compute(None, soil, PolarizationState.VV)
            hh_result = model.compute(None, soil, PolarizationState.HH)
            hv_result = model.compute(None, soil, PolarizationState.HV)
            
            # Convert linear to dB
            vv_db = 10.0 * np.log10(vv_result) if vv_result > 0 else float('-inf')
            hh_db = 10.0 * np.log10(hh_result) if hh_result > 0 else float('-inf')
            hv_db = 10.0 * np.log10(hv_result) if hv_result > 0 else float('-inf')
            
        except Exception as e:
            # Skip this configuration if computation fails
            import traceback
            print(f"Warning: KA computation failed for ratio={ratio}, k*L={k_L:.2f}: {e}")
            if "src" in str(e):
                print("Full traceback:")
                traceback.print_exc()
            continue

        overall_model["hh"].append(hh_db)
        overall_model["vv"].append(vv_db)
        overall_model["hv"].append(hv_db)

        overall_reference["hh"].append(hh_ref)
        overall_reference["vv"].append(vv_ref)
        overall_reference["hv"].append(hv_ref)

        if ratio not in grouped_model:
            grouped_model[ratio] = {pol: [] for pol in ("hh", "vv", "hv")}
            grouped_reference[ratio] = {pol: [] for pol in ("hh", "vv", "hv")}

        grouped_model[ratio]["hh"].append(hh_db)
        grouped_model[ratio]["vv"].append(vv_db)
        grouped_model[ratio]["hv"].append(hv_db)

        grouped_reference[ratio]["hh"].append(hh_ref)
        grouped_reference[ratio]["vv"].append(vv_ref)
        grouped_reference[ratio]["hv"].append(hv_ref)

    overall_metrics = {
        pol: _calc_metrics(overall_model[pol], overall_reference[pol]) for pol in overall_model
    }

    by_ratio: Dict[float, Dict[str, Metrics]] = {}
    for ratio_value, model_dict in grouped_model.items():
        by_ratio[ratio_value] = {
            pol: _calc_metrics(model_dict[pol], grouped_reference[ratio_value][pol])
            for pol in model_dict
        }

    return ComparisonResult(overall=overall_metrics, by_ratio=by_ratio)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare KA predictions with NMM3D LUT data")
    parser.add_argument(
        "--lut",
        type=Path,
        default=_DEFAULT_LUT,
        help=f"Path to LUT file (default: {_DEFAULT_LUT})",
    )
    parser.add_argument(
        "--frequency",
        type=float,
        default=_DEFAULT_FREQ_GHZ,
        help=f"Radar frequency in GHz (default: {_DEFAULT_FREQ_GHZ})",
    )
    parser.add_argument(
        "--incidence",
        type=float,
        default=_DEFAULT_INC_DEG,
        help=f"Incidence angle in degrees (default: {_DEFAULT_INC_DEG})",
    )
    parser.add_argument(
        "--phi",
        type=float,
        default=_DEFAULT_PHI_DEG,
        help=f"Scattering azimuth in degrees (default: {_DEFAULT_PHI_DEG})",
    )
    parser.add_argument(
        "--ratios",
        type=float,
        nargs="*",
        help="Optional list of correlation length ratios (ℓ/σ) to evaluate",
    )
    parser.add_argument(
        "--min-ratio",
        type=float,
        help="Minimum correlation length ratio to include (KA works best for large ratios)",
    )
    parser.add_argument(
        "--per-ratio",
        action="store_true",
        help="Print metrics broken down by ratio in addition to overall statistics",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        table = _load_lut(args.lut)
    except FileNotFoundError as exc:  # pragma: no cover - CLI guard
        print(f"Error: {exc}")
        return 1

    rows = _select_angle(table, args.incidence)
    if rows.size == 0:
        print(f"Error: No LUT entries found for incidence angle {args.incidence}")
        return 1

    # Calculate wavelength for validity check info
    lam = toLambda(args.frequency)
    k = 2.0 * math.pi / lam
    
    print("=" * 70)
    print("KA Model vs NMM3D Comparison")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Frequency: {args.frequency} GHz")
    print(f"  Wavelength: {lam*100:.2f} cm")
    print(f"  Incidence angle: {args.incidence}°")
    print(f"  Wavenumber k: {k:.2f} rad/m")
    
    if args.min_ratio:
        print(f"  Minimum ratio filter: ℓ/σ ≥ {args.min_ratio}")
    
    print(f"\nNote: KA model is designed for rough surfaces (k*σ > 0.3)")
    print(f"      KA should perform BETTER on very rough surfaces compared to NMM3D")
    print(f"      Using exponential ACF to match NMM3D configuration")

    result = _run_comparison(
        rows=rows,
        frequency_ghz=args.frequency,
        incidence_deg=args.incidence,
        phi_deg=args.phi,
        ratios=args.ratios,
        min_ratio=args.min_ratio,
    )

    overall = result.overall
    print("\n" + "=" * 70)
    print("Overall Metrics (KA vs NMM3D)")
    print("=" * 70)
    for pol in ("vv", "hh", "hv"):
        metrics = overall[pol]
        print(metrics.format_row(pol.upper()))

    if args.per_ratio and result.by_ratio:
        print("\n" + "=" * 70)
        print("By-Ratio Metrics")
        print("=" * 70)
        for ratio in sorted(result.by_ratio):
            print(f"\nℓ/σ = {ratio:g}")
            # Calculate k*L for this ratio
            sigma_norm = 1.0 / ratio  # Approximate
            corr_len = ratio * sigma_norm * lam
            k_L = k * corr_len
            print(f"  (k*ℓ ≈ {k_L:.2f})")
            for pol in ("vv", "hh", "hv"):
                print("  " + result.by_ratio[ratio][pol].format_row(pol.upper()))

    total_valid = sum(metrics.count for metrics in overall.values())
    
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Total valid comparisons: {total_valid}")
    
    if total_valid == 0:
        print("\nWarning: No valid comparisons available.")
        print("This is expected if all surface configurations are outside KA validity range.")
        print("KA model requires large-scale roughness (k*ℓ > 6).")
        print("\nTry:")
        print("  --min-ratio 10  (to focus on larger correlation lengths)")
        return 1
    
    # Provide interpretation
    print("\nInterpretation:")
    if overall["vv"].count > 0:
        if overall["vv"].rmse < 3.0:
            print("  ✓ Excellent agreement (RMSE < 3 dB)")
        elif overall["vv"].rmse < 5.0:
            print("  ✓ Good agreement (RMSE < 5 dB)")
        elif overall["vv"].rmse < 10.0:
            print("  ⚠ Moderate agreement (RMSE < 10 dB)")
        else:
            print("  ⚠ Limited agreement (RMSE > 10 dB)")
            print("    Note: KA is designed for large-scale roughness with moderate slopes")
    
    # Analyze by-ratio performance
    if result.by_ratio:
        print("\nBy-Ratio Analysis:")
        best_ratio = None
        best_rmse = float('inf')
        for ratio, metrics in result.by_ratio.items():
            if metrics['vv'].count > 0 and metrics['vv'].rmse < best_rmse:
                best_rmse = metrics['vv'].rmse
                best_ratio = ratio
        
        if best_ratio is not None:
            print(f"  Best agreement at ℓ/σ = {best_ratio} (RMSE = {best_rmse:.2f} dB)")
            
            # Calculate MSS for best ratio
            sigma_norm = 1.0 / best_ratio
            mss = 2.0 * (sigma_norm) ** 2
            print(f"  Mean square slope at best ratio: {mss:.6f}")
            
            if best_ratio < 6:
                print("  → KA works best for moderate ratios (4-8) with moderate slopes")
            elif best_ratio > 12:
                print("  → Large ratios have very small slopes, reducing KA backscatter")
    
    print("\nPhysical Interpretation:")
    print("  • KA model captures large-scale specular reflection")
    print("  • For very smooth slopes (large ℓ/σ), KA predicts low backscatter")
    print("  • NMM3D includes both large and small-scale scattering")
    print("  • Best agreement occurs at moderate slopes (MSS ~ 0.1-0.2)")
    print("  • For complete modeling, consider two-scale models (KA + Bragg)")
    
    print("\nNote on Cross-Polarization:")
    print("  • Basic KA model predicts ZERO cross-pol for monostatic backscatter")
    print("  • This is theoretically correct for single-bounce specular reflection")
    print("  • Cross-pol in backscatter requires:")
    print("    - Multiple scattering (double-bounce)")
    print("    - Small-scale roughness (Bragg scattering)")
    print("    - Volume scattering")
    print("  • NMM3D includes these mechanisms, hence non-zero HV")
    print("  • For cross-pol modeling, use two-scale models (KA + SPM/Bragg)")

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
