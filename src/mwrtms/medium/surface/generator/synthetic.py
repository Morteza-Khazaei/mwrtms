"""Synthetic self-affine surface generation."""

from __future__ import annotations

from typing import Tuple

import numpy as np

from ..base import Surface
from .base import SurfaceGenerator


class SyntheticSurfaceGenerator(SurfaceGenerator):
    """Generate synthetic self-affine rough surfaces using spectral synthesis."""

    def __init__(
        self,
        shape: Tuple[int, int],
        physical_size: Tuple[float, float],
        *,
        hurst: float = 0.8,
        rms: float = 1e-6,
        correlation_type: str = "exponential",
    ) -> None:
        super().__init__(shape, physical_size)

        if not (0.0 < hurst < 1.0):
            raise ValueError("hurst must lie in (0, 1)")
        if rms <= 0.0:
            raise ValueError("rms must be positive")

        self._hurst = float(hurst)
        self._rms = float(rms)
        self._correlation_type = correlation_type

    # Fluent configuration helpers -------------------------------------------------
    def setHurstExponent(self, hurst: float) -> "SyntheticSurfaceGenerator":  # noqa: N802
        self._hurst = float(hurst)
        return self

    def setRMSHeight(self, rms: float) -> "SyntheticSurfaceGenerator":  # noqa: N802
        if rms <= 0.0:
            raise ValueError("rms must be positive")
        self._rms = float(rms)
        return self

    def setCorrelationType(self, corr_type: str) -> "SyntheticSurfaceGenerator":  # noqa: N802
        self._correlation_type = corr_type
        return self

    # Internal generation ----------------------------------------------------------
    def _generate_impl(self) -> Surface:
        nx, ny = self._shape
        lx, ly = self._physical_size

        qx = np.fft.fftfreq(nx, d=lx / nx) * 2.0 * np.pi
        qy = np.fft.fftfreq(ny, d=ly / ny) * 2.0 * np.pi
        Qx, Qy = np.meshgrid(qx, qy, indexing="ij")
        Q = np.sqrt(Qx**2 + Qy**2)

        Q[0, 0] = np.finfo(float).eps

        spectral_exponent = -2.0 - 2.0 * self._hurst
        spectrum = np.power(Q, spectral_exponent)

        if self._correlation_type == "gaussian":
            cutoff = 2.0 * np.pi / min(lx, ly)
            spectrum *= np.exp(-(Q**2) / (2.0 * cutoff**2))
        elif self._correlation_type not in {"exponential", "gaussian"}:
            # Allow custom handling by raising early
            raise ValueError(f"unsupported correlation_type '{self._correlation_type}'")

        spectrum[0, 0] = 0.0

        phase = np.random.uniform(0.0, 2.0 * np.pi, size=(nx, ny))
        amplitudes = np.sqrt(spectrum)
        field_fft = amplitudes * (np.cos(phase) + 1j * np.sin(phase))

        heights = np.fft.ifft2(field_fft).real
        std = np.std(heights)
        if std == 0.0:
            raise RuntimeError("unable to generate non-flat surface")
        heights *= self._rms / std

        metadata = {
            "generator": "synthetic",
            "hurst": self._hurst,
            "correlation_type": self._correlation_type,
            "target_rms": self._rms,
        }

        return Surface(heights, self._physical_size, metadata)
