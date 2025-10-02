"""Loader for measured surface data."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple, Union

import numpy as np

from ..base import Surface
from .base import SurfaceGenerator


class MeasuredSurfaceLoader(SurfaceGenerator):
    """Generate a :class:`Surface` instance from measured height data."""

    def __init__(self, source: Union[np.ndarray, str, Path], physical_size: Tuple[float, float]) -> None:
        self._source = source
        heights = self._load_data()
        super().__init__(heights.shape, physical_size)
        self._heights = heights

    def _load_data(self) -> np.ndarray:
        if isinstance(self._source, np.ndarray):
            data = np.asarray(self._source, dtype=float)
            if data.ndim != 2:
                raise ValueError("measured surface array must be 2D")
            return data.copy()

        path = Path(self._source)
        if not path.exists():
            raise FileNotFoundError(path)

        suffix = path.suffix.lower()
        if suffix == ".npy":
            data = np.load(path)
        elif suffix in {".txt", ".csv"}:
            data = np.loadtxt(path)
        elif suffix == ".mat":
            try:
                from scipy.io import loadmat
            except ImportError as exc:  # pragma: no cover - SciPy should be available
                raise RuntimeError("scipy is required to load MAT files") from exc
            loaded = loadmat(path)
            numeric = [val for val in loaded.values() if isinstance(val, np.ndarray) and val.ndim == 2]
            if not numeric:
                raise ValueError(f"no 2D numeric arrays found in {path}")
            data = numeric[0]
        else:
            raise ValueError(f"unsupported file format '{suffix}'")

        data = np.asarray(data, dtype=float)
        if data.ndim != 2:
            raise ValueError("loaded data must be 2D")
        return data

    def _generate_impl(self) -> Surface:
        metadata = {
            "generator": "measured",
            "source": str(self._source) if not isinstance(self._source, np.ndarray) else "array",
        }
        return Surface(self._heights, self._physical_size, metadata)
