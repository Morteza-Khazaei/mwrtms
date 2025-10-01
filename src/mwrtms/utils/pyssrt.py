"""Utility helpers bridging mwRTMs objects with legacy pySSRT routines."""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np

from ..interface.correlation import ExponentialCorrelation, GaussianCorrelation, PowerLawCorrelation
from ..interface.roughness import SurfaceRoughness

__all__ = [
    "determine_acf_descriptor",
    "ensure_polarization_dict",
    "db_to_power",
]


def determine_acf_descriptor(
    roughness: SurfaceRoughness,
    override: str | None = None,
) -> tuple[str, float | None]:
    """Return the pySSRT autocorrelation identifier and optional exponent.

    Parameters
    ----------
    roughness:
        Roughness description supplying the correlation function and length.
    override:
        Optional manual label (``"gauss"``, ``"exp"``, ``"pow"``). When omitted
        the function infers the label from ``roughness.correlation_function``.

    Returns
    -------
    tuple
        ``(label, extra)`` suitable for the legacy routines where ``extra`` is
        the power-law exponent when required.
    """

    if override is not None:
        label = override.lower()
        extra = None
        if label == "pow":
            corr = roughness.correlation_function
            if isinstance(corr, PowerLawCorrelation):
                extra = corr.power
        return label, extra

    corr = roughness.correlation_function
    if isinstance(corr, GaussianCorrelation):
        return "gauss", None
    if isinstance(corr, ExponentialCorrelation):
        return "exp", None
    if isinstance(corr, PowerLawCorrelation):
        return "pow", corr.power
    # Default to exponential correlation for custom implementations.
    return "exp", None


def ensure_polarization_dict(values: Dict[str, float] | Iterable[Tuple[str, float]]) -> Dict[str, float]:
    """Return a lowercase-polarisation dictionary with missing keys filled.

    The helper standardises keys expected by mwRTMs (``"vv"``, ``"hh"``,
    ``"hv"``, ``"vh"``) and guarantees float values. Missing entries default to
    zero.
    """

    if isinstance(values, dict):
        raw = values.items()
    else:
        raw = values
    data: Dict[str, float] = {k.lower(): float(v) for k, v in raw}
    for key in ("vv", "hh", "hv", "vh"):
        data.setdefault(key, 0.0)
    return data


def db_to_power(value: float | np.ndarray) -> float | np.ndarray:
    """Convert backscatter expressed in dB to linear power."""

    return np.power(10.0, np.asarray(value, dtype=float) / 10.0)
