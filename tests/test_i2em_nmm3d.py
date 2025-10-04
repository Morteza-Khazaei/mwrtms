"""CLI harness comparing I2EM against the NMM3D LUT."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np

from mwrtms import (
    RadarConfigurationFactory,
    ElectromagneticWave,
    ScatteringGeometry,
    PolarizationState,
)
from mwrtms.core.constants import SPEED_OF_LIGHT
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface import I2EMModel

_DEFAULT_LUT = Path("data/NMM3D_LUT_NRCS_40degree.dat")
_DEFAULT_FREQ_GHZ = 5.405
_DEFAULT_INC_DEG = 40.0


@dataclass
class Metrics:
    count: int
    rmse: float
    mae: float
    bias: float
    corr: float

    def format_row(self, label: str) -> str:
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


def _calc_metrics(model: Iterable[float], reference: Iterable[float]) -> Metrics:
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
    overall: Dict[str, Metrics]
    by_ratio: Dict[float, Dict[str, Metrics]]


def _run_comparison(
    rows: np.ndarray,
    frequency_ghz: float,
    ratios: Sequence[float] | None,
) -> ComparisonResult:
    wavelength = SPEED_OF_LIGHT / (frequency_ghz * 1e9)

    wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
    geometry = ScatteringGeometry(theta_i_deg=_DEFAULT_INC_DEG)
    air = HomogeneousMedium(1.0 + 0.0j)

    overall_model = {pol: [] for pol in ("vv", "hh")}
    overall_ref = {pol: [] for pol in ("vv", "hh")}
    grouped_model: Dict[float, Dict[str, List[float]]] = {}
    grouped_ref: Dict[float, Dict[str, List[float]]] = {}

    for row in rows:
        _, ratio, eps_real, eps_imag, rms_norm, vv_ref, hh_ref, _ = row
        ratio = float(ratio)
        if ratios and ratio not in ratios:
            continue

        sigma = rms_norm * wavelength
        corr = ratio * sigma
        surface = build_surface_from_statistics(sigma, corr, correlation_type="exponential")
        soil = HomogeneousMedium(complex(eps_real, eps_imag))

        model = I2EMModel(wave, geometry, surface, correlation_type="exponential")
        result = model.run(air, soil)

        vv_db = 10.0 * np.log10(result[PolarizationState.VV] + 1e-30)
        hh_db = 10.0 * np.log10(result[PolarizationState.HH] + 1e-30)

        overall_model["vv"].append(vv_db)
        overall_model["hh"].append(hh_db)
        overall_ref["vv"].append(vv_ref)
        overall_ref["hh"].append(hh_ref)

        grouped_model.setdefault(ratio, {"vv": [], "hh": []})
        grouped_ref.setdefault(ratio, {"vv": [], "hh": []})
        grouped_model[ratio]["vv"].append(vv_db)
        grouped_model[ratio]["hh"].append(hh_db)
        grouped_ref[ratio]["vv"].append(vv_ref)
        grouped_ref[ratio]["hh"].append(hh_ref)

    overall_metrics = {pol: _calc_metrics(overall_model[pol], overall_ref[pol]) for pol in overall_model}
    by_ratio = {
        ratio: {pol: _calc_metrics(grouped_model[ratio][pol], grouped_ref[ratio][pol]) for pol in grouped_model[ratio]}
        for ratio in grouped_model
    }

    return ComparisonResult(overall=overall_metrics, by_ratio=by_ratio)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare I2EM predictions with NMM3D LUT data")
    parser.add_argument("--lut", type=Path, default=_DEFAULT_LUT, help="Path to LUT file")
    parser.add_argument("--frequency", type=float, default=_DEFAULT_FREQ_GHZ, help="Radar frequency in GHz")
    parser.add_argument(
        "--ratios",
        type=float,
        nargs="*",
        help="Optional list of correlation ratios to evaluate",
    )
    parser.add_argument("--per-ratio", action="store_true", help="Print metrics broken down by ratio")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    table = _load_lut(args.lut)
    rows = _select_angle(table, _DEFAULT_INC_DEG)
    if rows.size == 0:
        parser.error(f"No LUT entries found at {_DEFAULT_INC_DEG} degrees")

    result = _run_comparison(rows, args.frequency, args.ratios)

    print("I2EM vs NMM3D (overall metrics)")
    for pol in ("vv", "hh"):
        print(result.overall[pol].format_row(pol.upper()))

    if args.per_ratio:
        print("\nBy-ratio metrics")
        for ratio in sorted(result.by_ratio):
            print(f"\nℓ/σ = {ratio:g}")
            for pol in ("vv", "hh"):
                print(result.by_ratio[ratio][pol].format_row(pol.upper()))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
