"""Surface analysis utilities."""

from __future__ import annotations

from typing import Dict

import numpy as np

from .base import Surface


class SurfaceAnalyzer:
    """High level helpers to derive parameters from :class:`Surface`."""

    @staticmethod
    def analyze(surface: Surface) -> Dict[str, float]:
        """Return a dictionary of key surface statistics."""

        results: Dict[str, float] = {
            "rms_height": surface.rms_height(),
            "mean_height": surface.mean_height(),
            "peak_to_valley": surface.peak_to_valley(),
            "rms_slope": surface.rms_slope(),
            "is_isotropic": surface.is_isotropic(),
            "anisotropy_ratio": surface.anisotropy_ratio(),
            "correlation_length": surface.correlation_length(method="acf"),
            "correlation_length_spectral": surface.correlation_length(method="spectral"),
        }

        if not results["is_isotropic"]:
            ell_x, ell_y = surface.correlation_lengths_xy()
            results["correlation_length_x"] = ell_x
            results["correlation_length_y"] = ell_y

        return results

    @staticmethod
    def estimate_correlation_type(surface: Surface) -> str:
        """Estimate the best fitting correlation function type."""

        metadata = surface.metadata
        corr_meta = metadata.get("correlation_type")
        if isinstance(corr_meta, str):
            return corr_meta

        try:

            from scipy.optimize import curve_fit
        except ImportError:  # pragma: no cover - SciPy should be present
            return "exponential"

        lags, acf = surface.radial_autocorrelation()
        if lags.size < 5:
            return "exponential"

        # Ignore the first entry (zero lag) to avoid singularities in fitting.
        lags_fit = lags[1: min(len(lags), 50)]
        acf_fit = acf[1: min(len(acf), 50)]

        if lags_fit.size < 3:
            return "exponential"

        def exp_model(r, ell):
            return np.exp(-r / np.maximum(ell, 1e-12))

        def gauss_model(r, ell):
            return np.exp(-((r / np.maximum(ell, 1e-12)) ** 2))

        threshold = np.exp(-1.0)
        below = np.where(acf_fit < threshold)[0]
        if below.size:
            initial = float(max(lags_fit[below[0]], 1e-12))
        else:
            initial = float(max(lags_fit[len(lags_fit) // 2], 1e-12))
        p0 = (initial,)

        def r_squared(model) -> float:
            try:
                popt, _ = curve_fit(model, lags_fit, acf_fit, p0=p0, maxfev=5000)
            except Exception:  # pragma: no cover - robust fallback
                return -np.inf
            predictions = model(lags_fit, *popt)
            residuals = acf_fit - predictions
            ss_res = float(np.sum(residuals**2))
            ss_tot = float(np.sum((acf_fit - np.mean(acf_fit)) ** 2))
            if ss_tot == 0.0:
                return -np.inf
            return 1.0 - ss_res / ss_tot

        r2_exp = r_squared(exp_model)
        r2_gauss = r_squared(gauss_model)

        if r2_gauss - r2_exp > 0.05:
            return "gaussian"
        return "exponential"

    @staticmethod
    def correlation_model(surface: Surface):
        """Return (correlation_function, lengths) tuple for ``surface``."""

        corr_type = SurfaceAnalyzer.estimate_correlation_type(surface)
        from .correlation import (
            Exponential,
            Gaussian,
            PowerLaw,
        )

        if corr_type == "gaussian":
            corr_func = Gaussian()
        elif corr_type == "powerlaw":
            corr_func = PowerLaw()
        else:
            corr_func = Exponential()

        if surface.is_isotropic():
            lengths = surface.correlation_length()
        else:
            lengths = surface.correlation_lengths_xy()

        return corr_func, lengths
