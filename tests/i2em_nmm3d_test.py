"""CLI test harness comparing I2EM against the NMM3D LUT.

This script provides a command-line interface to perform numerical checks
similar to those in the `test_i2em.ipynb` notebook. It is designed for
lightweight, automation-friendly validation. The script loads the 40° incidence
NMM3D backscatter look-up table, evaluates the I2EM model for each corresponding
surface configuration, and reports goodness-of-fit metrics (RMSE, MAE, bias,
and Pearson correlation) for the HH, VV, and HV channels.

Usage
-----
Run from the repository root:

    python3 tests/i2em_nmm3d_test.py

Optional command-line arguments allow for filtering by ratio, specifying a
different LUT file, or changing the radar frequency. The script will exit with
a non-zero status if the LUT is not found or if no valid comparisons can be
made after applying filters.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

from mwrtms.core import ElectromagneticWave, PolarizationState, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface import I2EMModel

# Default configuration mirrors the notebook and other tests
_DEFAULT_LUT = Path(__file__).resolve().parents[1] / "data" / "NMM3D_LUT_NRCS_40degree.dat"
_DEFAULT_FREQ_GHZ = 5.405
_DEFAULT_INC_DEG = 40.0


@dataclass
class Metrics:
    """Container for descriptive statistics between I2EM and NMM3D."""

    count: int
    rmse: float
    mae: float
    bias: float
    corr: float

    def format_row(self, label: str) -> str:
        """Format metrics into a standardized string for printing."""
        return (
            f"{label:<6s} n={self.count:3d}  "
            f"RMSE={self.rmse:5.2f} dB  "
            f"MAE={self.mae:5.2f} dB  "
            f"Bias={self.bias:+5.2f} dB  "
            f"Corr={self.corr:6.3f}"
        )


def _to_db(value: float) -> float:
    """Safely convert a linear value to decibels."""
    tiny = np.finfo(float).tiny
    return 10.0 * np.log10(max(value, tiny))


def _load_lut(path: Path) -> np.ndarray:
    """Load the NMM3D LUT from the specified path."""
    if not path.exists():
        raise FileNotFoundError(f"LUT not found at {path}")
    return np.loadtxt(path)


def _select_angle(table: np.ndarray, incidence_deg: float, tol: float = 1e-3) -> np.ndarray:
    """Filter LUT rows by a specific incidence angle."""
    mask = np.abs(table[:, 0] - incidence_deg) <= tol
    return table[mask]


def _calc_metrics(model: Iterable[float], reference: Iterable[float]) -> Metrics:
    """Calculate performance metrics between model and reference data."""
    model_arr = np.asarray(list(model), dtype=float)
    ref_arr = np.asarray(list(reference), dtype=float)
    mask = np.isfinite(model_arr) & np.isfinite(ref_arr)
    if not np.any(mask):
        return Metrics(0, float("nan"), float("nan"), float("nan"), float("nan"))

    diff = model_arr[mask] - ref_arr[mask]
    rmse = float(np.sqrt(np.mean(diff**2)))
    mae = float(np.mean(np.abs(diff)))
    bias = float(np.mean(diff))
    corr = float(np.corrcoef(ref_arr[mask], model_arr[mask])[0, 1]) if mask.sum() > 1 else float("nan")
    return Metrics(int(mask.sum()), rmse, mae, bias, corr)


@dataclass
class ComparisonResult:
    """Holds the results of the comparison, both overall and grouped by ratio."""

    overall: Dict[str, Metrics]
    by_ratio: Dict[float, Dict[str, Metrics]]


def _run_comparison(
    rows: np.ndarray,
    frequency_ghz: float,
    incidence_deg: float,
    ratios: Sequence[float] | None,
) -> ComparisonResult:
    """Execute the I2EM vs. NMM3D comparison."""
    wave = ElectromagneticWave(frequency_ghz * 1e9)
    geometry = ScatteringGeometry(theta_i_deg=incidence_deg)
    air = HomogeneousMedium(1.0 + 0.0j)
    lambda_m = wave.wavelength

    results: Dict[float, Dict[str, List[float]]] = {}

    for row in rows:
        _, ratio, eps_r, eps_i, rms_norm, vv_ref, hh_ref, hv_ref = row
        ratio = float(ratio)
        if ratios and ratio not in ratios:
            continue

        sigma = rms_norm * lambda_m
        surface = build_surface_from_statistics(sigma, ratio * sigma, correlation_type="exponential")
        soil = HomogeneousMedium(complex(eps_r, eps_i))
        model = I2EMModel(wave, geometry, surface)

        results.setdefault(ratio, {"vv_i2em": [], "hh_i2em": [], "hv_i2em": [], "vv_ref": [], "hh_ref": [], "hv_ref": []})
        results[ratio]["vv_i2em"].append(_to_db(model.compute(air, soil, PolarizationState.VV)))
        results[ratio]["hh_i2em"].append(_to_db(model.compute(air, soil, PolarizationState.HH)))
        if np.isfinite(hv_ref):
            results[ratio]["hv_i2em"].append(_to_db(model.compute(air, soil, PolarizationState.HV)))
        else:
            results[ratio]["hv_i2em"].append(float("-inf"))
        results[ratio]["vv_ref"].append(vv_ref)
        results[ratio]["hh_ref"].append(hh_ref)
        results[ratio]["hv_ref"].append(hv_ref)

    overall_metrics = {
        pol: _calc_metrics(
            np.concatenate([r[f"{pol}_i2em"] for r in results.values()]),
            np.concatenate([r[f"{pol}_ref"] for r in results.values()]),
        )
        for pol in ("vv", "hh", "hv")
    }

    by_ratio_metrics = {
        ratio: {pol: _calc_metrics(res[f"{pol}_i2em"], res[f"{pol}_ref"]) for pol in ("vv", "hh", "hv")}
        for ratio, res in results.items()
    }

    return ComparisonResult(overall=overall_metrics, by_ratio=by_ratio_metrics)


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(description="Compare I2EM predictions with NMM3D LUT data")
    parser.add_argument("--lut", type=Path, default=_DEFAULT_LUT, help=f"Path to LUT file (default: {_DEFAULT_LUT})")
    parser.add_argument("--frequency", type=float, default=_DEFAULT_FREQ_GHZ, help=f"Radar frequency in GHz (default: {_DEFAULT_FREQ_GHZ})")
    parser.add_argument("--ratios", type=float, nargs="*", help="Optional list of correlation ratios (ℓ/σ) to evaluate")
    parser.add_argument("--per-ratio", action="store_true", help="Print metrics broken down by ratio")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Main entry point for the script."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        table = _load_lut(args.lut)
    except FileNotFoundError as exc:
        parser.error(str(exc))

    rows = _select_angle(table, _DEFAULT_INC_DEG)
    if rows.size == 0:
        parser.error(f"No LUT entries found at {_DEFAULT_INC_DEG} degrees")

    result = _run_comparison(rows, args.frequency, _DEFAULT_INC_DEG, args.ratios)

    print("I2EM vs NMM3D (overall metrics)")
    for pol in ("vv", "hh", "hv"):
        print(result.overall[pol].format_row(pol.upper()))

    if args.per_ratio:
        print("\nBy-ratio metrics")
        for ratio in sorted(result.by_ratio):
            print(f"\nℓ/σ = {ratio:g}")
            for pol in ("vv", "hh", "hv"):
                print(result.by_ratio[ratio][pol].format_row(pol.upper()))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
