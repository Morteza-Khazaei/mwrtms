"""Base utilities for IEM-family scattering models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Union

import numpy as np

from ..base import SurfaceScattering
from ....core import PolarizationState
from ....core.geometry import ScatteringGeometry
from ....core.polarization import normalize_polarization
from ....core.radar_modes import RadarConfiguration
from ....core.wave import ElectromagneticWave
from ....medium import Medium
from ....result.scattering import ScatteringResult
from .correlation import (
    auto_spectral_order,
    single_scale_spectrum_weights,
)

PolarizationInput = Union[PolarizationState, str, Iterable[Union[PolarizationState, str]], None]


@dataclass
class SurfaceRoughnessParameters:
    sigma_m: float
    correlation_length_m: float
    ks: float
    kl: float
    rms_slope: float
    spectral_weights: np.ndarray


class IEMBase(SurfaceScattering):
    """Shared functionality for IEM-family surface scattering models."""

    MODEL_NAME = "IEM"

    def __init__(
        self,
        wave: ElectromagneticWave,
        geometry: ScatteringGeometry,
        surface,
        *,
        correlation_type: str = "exponential",
        correlation_length_m: Optional[float] = None,
        spectral_terms: Optional[int] = None,
        power_exponent: float = 1.5,
        auto_terms: bool = True,
    ) -> None:
        super().__init__(wave, geometry, surface=surface)
        self._correlation_type = correlation_type
        self._correlation_length_m = correlation_length_m
        self._spectral_terms = spectral_terms
        self._power_exponent = power_exponent
        self._auto_terms = auto_terms

    def _surface_parameters(self) -> SurfaceRoughnessParameters:
        sigma_m = self._surface.rms_height()
        corr_m = self._correlation_length_m or self._surface.correlation_length()
        k = self._wave.wavenumber
        ks = k * sigma_m
        kl = k * corr_m

        if self._auto_terms:
            terms = auto_spectral_order(ks, np.cos(self._geometry.theta_i_rad), np.cos(self._geometry.theta_s_rad))
        else:
            if self._spectral_terms is None:
                raise ValueError("spectral_terms must be provided when auto_terms is False")
            terms = self._spectral_terms

        wvnb = self._compute_roughness_wavenumber()
        weights, slope = single_scale_spectrum_weights(
            wvnb=wvnb,
            sigma=sigma_m,
            corr_length=corr_m,
            n_terms=terms,
            correlation=self._correlation_type,
            power_exponent=self._power_exponent,
        )
        rms = slope
        return SurfaceRoughnessParameters(
            sigma_m=sigma_m,
            correlation_length_m=corr_m,
            ks=ks,
            kl=kl,
            rms_slope=rms,
            spectral_weights=weights,
        )

    def _compute_roughness_wavenumber(self) -> float:
        k = self._wave.wavenumber
        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_s_rad
        phi_s = self._geometry.phi_s_rad
        return k * np.sqrt(
            (np.sin(theta_s) * np.cos(phi_s) - np.sin(theta_i)) ** 2
            + (np.sin(theta_s) * np.sin(phi_s)) ** 2
        )

    def run(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarizations: PolarizationInput = None,
        radar_config: Optional[RadarConfiguration] = None,
    ) -> ScatteringResult:
        params = self._surface_parameters()
        pol_tuple = normalize_polarization(polarizations)
        data = {
            pol: self._compute_channel(medium_above, medium_below, pol, params)
            for pol in pol_tuple
        }
        return ScatteringResult(
            data,
            self.MODEL_NAME,
            self._wave,
            self._geometry,
            radar_config=radar_config,
        )

    def compute(self, medium_above, medium_below, polarization) -> float:
        params = self._surface_parameters()
        pol = normalize_polarization(polarization)[0]
        return self._compute_channel(medium_above, medium_below, pol, params)

    def _compute_channel(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
        params: SurfaceRoughnessParameters,
    ) -> float:
        raise NotImplementedError("Subclasses must implement _compute_channel")

    # SurfaceScattering abstract hooks are unused for IEM family but
    # must be defined to satisfy the interface. They simply delegate to
    # the more general `run` implementation.

    def _compute_kirchhoff(self, R_h, R_v, polarization):  # pragma: no cover - not used
        return 0.0

    def _compute_complementary(self, R_h, R_v, polarization):  # pragma: no cover - not used
        return 0.0
