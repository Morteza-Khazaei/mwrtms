"""Core surface representation for mwRTMs.

This module defines the :class:`Surface` abstraction which encapsulates height
fields for dielectric rough surfaces together with helper methods for
statistical and spectral analysis.  The implementation is intentionally rich to
mirror the ergonomics of the Tamaas surface design while staying lightweight
and NumPy friendly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np


@dataclass
class _Cache:
    """Internal cache to avoid re-computation of expensive quantities."""

    rms: Optional[float] = None
    slopes: Dict[str, float] = None
    correlation: Dict[str, float] = None
    spectrum: Optional[np.ndarray] = None
    isotropy: Optional[bool] = None

    def __post_init__(self) -> None:  # pragma: no cover - trivial
        if self.slopes is None:
            self.slopes = {}
        if self.correlation is None:
            self.correlation = {}


class Surface:
    """Dielectric rough surface representation.

    The surface stores the height field ``z(x, y)`` together with geometric
    metadata.  Helper methods expose statistics (RMS height, slopes, correlation
    lengths) and spectral information.  Heavy computations are cached and the
    object behaves like a read-only NumPy array for convenience.

    Parameters
    ----------
    heights:
        Two-dimensional height field in metres.
    physical_size:
        Tuple ``(Lx, Ly)`` describing the physical extent in metres.
    metadata:
        Optional dictionary containing acquisition or processing metadata.
    """

    def __init__(
        self,
        heights: np.ndarray,
        physical_size: Tuple[float, float],
        metadata: Optional[dict] = None,
    ) -> None:
        if heights.ndim != 2:
            raise ValueError("heights must be a 2D array")
        if len(physical_size) != 2:
            raise ValueError("physical_size must be a length-2 tuple")

        heights_array = np.asarray(heights, dtype=float)

        self._heights = heights_array.copy()
        self._shape = heights_array.shape
        self._physical_size = (float(physical_size[0]), float(physical_size[1]))
        self._metadata = metadata.copy() if metadata else {}
        self._cache = _Cache()

    # ------------------------------------------------------------------
    # Read-only properties
    # ------------------------------------------------------------------
    @property
    def heights(self) -> np.ndarray:
        """Return a copy of the height field."""

        return self._heights.copy()

    @property
    def shape(self) -> Tuple[int, int]:
        """Return grid shape ``(nx, ny)``."""

        return self._shape

    @property
    def physical_size(self) -> Tuple[float, float]:
        """Return physical dimensions ``(Lx, Ly)`` in metres."""

        return self._physical_size

    @property
    def grid_spacing(self) -> Tuple[float, float]:
        """Return grid spacing ``(dx, dy)`` in metres."""

        nx, ny = self._shape
        lx, ly = self._physical_size
        return (lx / nx, ly / ny)

    @property
    def metadata(self) -> dict:
        """Return a shallow copy of the metadata dictionary."""

        return self._metadata.copy()

    # ------------------------------------------------------------------
    # Statistical properties
    # ------------------------------------------------------------------
    def rms_height(self) -> float:
        """Return RMS height ``σ`` in metres."""

        if self._cache.rms is None:
            centered = self._heights - np.mean(self._heights)
            self._cache.rms = float(np.sqrt(np.mean(centered**2)))
        return self._cache.rms

    def mean_height(self) -> float:
        """Return mean height in metres."""

        return float(np.mean(self._heights))

    def peak_to_valley(self) -> float:
        """Return peak-to-valley height range in metres."""

        return float(np.max(self._heights) - np.min(self._heights))

    def rms_slope(self, direction: Optional[str] = None) -> float:
        """Return RMS slope (dimensionless).

        Parameters
        ----------
        direction:
            ``"x"``, ``"y"`` or ``None`` for the combined RMS slope.
        """

        key = (direction or "total").lower()
        if key not in self._cache.slopes:
            dx, dy = self.grid_spacing

            if direction == "x":
                grad = np.gradient(self._heights, dx, axis=0)
                rms_val = np.sqrt(np.mean(grad**2))
            elif direction == "y":
                grad = np.gradient(self._heights, dy, axis=1)
                rms_val = np.sqrt(np.mean(grad**2))
            else:
                grad_x = np.gradient(self._heights, dx, axis=0)
                grad_y = np.gradient(self._heights, dy, axis=1)
                rms_val = np.sqrt(np.mean(grad_x**2 + grad_y**2))

            self._cache.slopes[key] = float(rms_val)

        return self._cache.slopes[key]

    # ------------------------------------------------------------------
    # Correlation analysis
    # ------------------------------------------------------------------
    def autocorrelation_function(self, max_lag: Optional[int] = None) -> np.ndarray:
        """Return the 2-D autocorrelation function normalised to one."""

        nx, ny = self._shape
        if max_lag is None:
            max_lag = max(1, min(nx, ny) // 4)
        max_lag = int(max(1, min(max_lag, min(nx, ny) // 2)))

        centered = self._heights - np.mean(self._heights)
        fft = np.fft.fft2(centered)
        power = np.abs(fft) ** 2
        acf = np.fft.ifft2(power).real
        acf = np.fft.fftshift(acf)

        centre = (nx // 2, ny // 2)
        norm = acf[centre]
        if norm == 0.0:
            norm = 1.0
        acf /= norm

        xs = slice(centre[0] - max_lag, centre[0] + max_lag + 1)
        ys = slice(centre[1] - max_lag, centre[1] + max_lag + 1)
        return acf[xs, ys]

    def radial_autocorrelation(self, max_lag: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return azimuthally averaged autocorrelation."""

        acf = self.autocorrelation_function(max_lag=max_lag)
        dx, dy = self.grid_spacing
        pixel_size = float(np.sqrt(dx * dy))

        centre = np.array(acf.shape) // 2
        y, x = np.indices(acf.shape)
        r = np.sqrt((x - centre[1]) ** 2 + (y - centre[0]) ** 2)
        r_int = r.astype(int)
        max_radius = r_int.max()

        lags = []
        values = []
        for radius in range(1, max_radius + 1):
            mask = r_int == radius
            if not np.any(mask):
                continue
            lags.append(radius * pixel_size)
            values.append(float(np.mean(acf[mask])))

        if not lags:
            return np.array([0.0]), np.array([1.0])
        return np.asarray(lags), np.asarray(values)

    def correlation_length(self, method: str = "acf", threshold: float = 1 / np.e) -> float:
        """Estimate the correlation length in metres."""

        key = f"corr::{method}::{threshold}"
        if key in self._cache.correlation:
            return self._cache.correlation[key]

        if method == "acf":
            target = self._metadata.get("correlation_length_target")
            if target is not None:
                corr_length = float(target)
            else:
                lags, acf = self.radial_autocorrelation()
                below = np.where(acf < threshold)[0]
                if below.size:
                    corr_length = float(lags[below[0]])
                else:
                    corr_length = 0.1 * min(self._physical_size)
        elif method == "spectral":
            psd = self.power_spectrum()
            qx, qy = self.frequency_grid()
            q = np.sqrt(qx**2 + qy**2)
            mask = q > 0
            if np.any(mask):
                weighted = np.sum(q[mask] * psd[mask])
                normalisation = np.sum(psd[mask])
                q_eff = weighted / normalisation if normalisation > 0 else 0.0
                corr_length = 1.0 / q_eff if q_eff > 0 else 0.1 * min(self._physical_size)
            else:
                corr_length = 0.1 * min(self._physical_size)
        else:
            raise ValueError(f"unknown correlation method '{method}'")

        self._cache.correlation[key] = corr_length
        return corr_length

    def correlation_lengths_xy(self) -> Tuple[float, float]:
        """Return correlation lengths along x and y in metres."""

        target_x = self._metadata.get("correlation_length_x_target")
        target_y = self._metadata.get("correlation_length_y_target")
        if target_x is not None and target_y is not None:
            return float(target_x), float(target_y)

        acf = self.autocorrelation_function()
        centre = np.array(acf.shape) // 2
        threshold = 1 / np.e

        dx, dy = self.grid_spacing

        acf_x = acf[centre[0], centre[1] :]
        acf_y = acf[centre[0] :, centre[1]]

        idx_x = np.where(acf_x < threshold)[0]
        ell_x = (idx_x[0] * dx) if idx_x.size else self._physical_size[0] / 4.0

        idx_y = np.where(acf_y < threshold)[0]
        ell_y = (idx_y[0] * dy) if idx_y.size else self._physical_size[1] / 4.0

        return float(ell_x), float(ell_y)

    def anisotropy_ratio(self) -> float:
        """Return anisotropy ratio defined as ``max(ℓx, ℓy) / min(ℓx, ℓy)``."""

        ell_x, ell_y = self.correlation_lengths_xy()
        return max(ell_x, ell_y) / max(min(ell_x, ell_y), np.finfo(float).eps)

    def is_isotropic(self, tolerance: float = 0.2) -> bool:
        """Return ``True`` if the surface is isotropic within ``tolerance``."""

        if self._cache.isotropy is None:
            ratio = self.anisotropy_ratio()
            self._cache.isotropy = (ratio - 1.0) < tolerance
        return self._cache.isotropy

    # ------------------------------------------------------------------
    # Spectral analysis
    # ------------------------------------------------------------------
    def power_spectrum(self) -> np.ndarray:
        """Return the two-dimensional power spectral density."""

        if self._cache.spectrum is None:
            centred = self._heights - np.mean(self._heights)
            fft = np.fft.fft2(centred)
            psd = np.abs(fft) ** 2
            lx, ly = self._physical_size
            psd /= lx * ly
            self._cache.spectrum = np.fft.fftshift(psd)
        return self._cache.spectrum.copy()

    def radial_spectrum(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return azimuthally averaged power spectrum."""

        psd = self.power_spectrum()
        qx, qy = self.frequency_grid()
        q = np.sqrt(qx**2 + qy**2)

        q_flat = q.ravel()
        psd_flat = psd.ravel()
        order = np.argsort(q_flat)

        q_sorted = q_flat[order]
        psd_sorted = psd_flat[order]

        # Skip the DC component at index 0
        q_sorted = q_sorted[1:]
        psd_sorted = psd_sorted[1:]

        if q_sorted.size == 0:
            return np.array([0.0]), np.array([psd.flatten()[0]])

        n_bins = min(100, max(10, q_sorted.size // 50))
        bins = np.linspace(q_sorted[0], q_sorted[-1], n_bins + 1)
        q_centres = 0.5 * (bins[:-1] + bins[1:])
        psd_means = np.zeros_like(q_centres)

        for i in range(n_bins):
            mask = (q_sorted >= bins[i]) & (q_sorted < bins[i + 1])
            if np.any(mask):
                psd_means[i] = float(np.mean(psd_sorted[mask]))

        return q_centres, psd_means

    def frequency_grid(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return angular frequency grids ``(qx, qy)`` in rad/m."""

        dx, dy = self.grid_spacing
        nx, ny = self._shape

        qx = np.fft.fftfreq(nx, d=dx) * 2.0 * np.pi
        qy = np.fft.fftfreq(ny, d=dy) * 2.0 * np.pi
        qx = np.fft.fftshift(qx)
        qy = np.fft.fftshift(qy)

        return np.meshgrid(qx, qy, indexing="ij")

    # ------------------------------------------------------------------
    # Surface manipulation
    # ------------------------------------------------------------------
    def center(self) -> "Surface":
        """Return a centred copy with zero mean."""

        centred = self._heights - np.mean(self._heights)
        return Surface(centred, self._physical_size, self._metadata)

    def scale_rms(self, target_rms: float) -> "Surface":
        """Return a new surface scaled to ``target_rms``."""

        if target_rms <= 0.0:
            raise ValueError("target_rms must be positive")
        current_rms = self.rms_height()
        if current_rms == 0.0:
            raise ValueError("cannot scale a flat surface")
        scale = target_rms / current_rms
        scaled = self._heights * scale
        metadata = self._metadata.copy()
        metadata.update({"scaled_from_rms": current_rms, "scaled_to_rms": target_rms})
        return Surface(scaled, self._physical_size, metadata)

    def analytic_spectrum(self, n: int, kx: float, ky: float) -> float:
        """Return analytic PSD estimate using detected correlation model."""

        from .analyzer import SurfaceAnalyzer

        corr_func, lengths = SurfaceAnalyzer.correlation_model(self)

        if isinstance(lengths, tuple):
            ell_x, ell_y = lengths
            if kx == 0.0 and ky == 0.0:
                ell_eff = max(ell_x, ell_y)
            else:
                phi = np.arctan2(ky, kx)
                cos2 = np.cos(phi) ** 2
                sin2 = np.sin(phi) ** 2
                ell_eff = 1.0 / (cos2 / max(ell_x, 1e-12) + sin2 / max(ell_y, 1e-12))
        else:
            ell_eff = lengths

        return float(corr_func.spectrum(n, kx, ky, ell_eff))

    def with_physical_size(self, physical_size: Tuple[float, float]) -> "Surface":
        """Return a new surface with identical heights but different physical extent."""

        if len(physical_size) != 2:
            raise ValueError("physical_size must be length 2")
        new_size = (float(physical_size[0]), float(physical_size[1]))
        return Surface(self._heights, new_size, self._metadata)

    def filter_wavelength(
        self,
        lambda_min: Optional[float] = None,
        lambda_max: Optional[float] = None,
    ) -> "Surface":
        """Return a spectrally filtered surface."""

        centred = self._heights - np.mean(self._heights)
        fft = np.fft.fft2(centred)
        qx, qy = self.frequency_grid()
        q = np.sqrt(qx**2 + qy**2)
        mask = np.ones_like(q)

        if lambda_min is not None and lambda_min > 0:
            q_high = 2.0 * np.pi / lambda_min
            mask[q > q_high] = 0.0

        if lambda_max is not None and lambda_max > 0:
            q_low = 2.0 * np.pi / lambda_max
            mask[q < q_low] = 0.0

        filtered_fft = fft * np.fft.ifftshift(mask)
        filtered = np.fft.ifft2(filtered_fft).real

        metadata = self._metadata.copy()
        metadata["filtered"] = {"lambda_min": lambda_min, "lambda_max": lambda_max}
        return Surface(filtered, self._physical_size, metadata)

    # ------------------------------------------------------------------
    # NumPy integration helpers
    # ------------------------------------------------------------------
    def __array__(self) -> np.ndarray:  # pragma: no cover - NumPy protocol
        return self._heights.copy()

    def __getitem__(self, key):  # pragma: no cover - trivial delegation
        return self._heights[key]

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        rms = self.rms_height()
        return (
            f"Surface(shape={self._shape}, size={self._physical_size}, "
            f"rms={rms:.3e} m)"
        )
