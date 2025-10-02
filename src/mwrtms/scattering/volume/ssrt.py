"""Single-Scattering Radiative Transfer (SSRT) vegetation model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...factory import register_model
from .base import VolumeScattering

__all__ = ["CanopyProperties", "SSRTModel"]


@dataclass
class CanopyProperties:
    """Physical description of the vegetation canopy."""

    lai: float
    vegetation_water_content: float
    height: float
    leaf_length: float = 0.05
    leaf_width: float = 0.02

    def __post_init__(self) -> None:
        if self.lai < 0.0:
            raise ValueError("LAI must be non-negative")
        if self.height <= 0.0:
            raise ValueError("Height must be positive")


@register_model("ssrt")
class SSRTModel(VolumeScattering):
    """Single-Scattering Radiative Transfer model with Rayleigh layer."""

    MODEL_NAME = "SSRT"

    def __init__(self, wave, geometry, canopy_properties: CanopyProperties, soil_backscatter: dict[str, float]) -> None:
        super().__init__(wave, geometry)
        self._canopy = canopy_properties
        self._soil_sigma0 = {}
        for key, value in soil_backscatter.items():
            if isinstance(key, str):
                self._soil_sigma0[key.lower()] = float(value)
            else:
                self._soil_sigma0[getattr(key, "value", str(key)).lower()] = float(value)

    def compute(self, medium_above, medium_below, polarization) -> float:
        sigma_vol = self._compute_volume_scattering(medium_below, polarization)
        sigma_gc = 0.0
        sigma_cg = 0.0
        tau_sq = self._compute_two_way_attenuation(medium_below, polarization)
        soil_sigma = self._soil_sigma0.get(polarization.value, 0.0)
        sigma_soil_att = tau_sq * soil_sigma
        return sigma_vol + sigma_gc + sigma_cg + sigma_soil_att

    def _compute_volume_scattering(self, veg_medium, polarization) -> float:
        k = self._wave.wavenumber
        theta = self._geometry.theta_i_rad
        eps_v = veg_medium.permittivity(self._wave.frequency_hz)
        a = self._canopy.leaf_length / 2.0
        ka = k * a
        if ka >= 0.3:
            import warnings
            warnings.warn(f"Rayleigh approximation limit reached: ka={ka:.3f}")
        volume = (4.0 / 3.0) * np.pi * a**3
        K = (eps_v - 1.0) / (eps_v + 2.0)
        sigma_leaf = (k**4 * volume**2 / np.pi) * np.abs(K) ** 2
        LAI = self._canopy.lai
        height = self._canopy.height
        n_leaf = LAI / height if height > 0 else 0.0
        kappa_e = 0.5 * LAI
        z_layers = np.linspace(0.0, height, 20)
        dz = height / 20 if height > 0 else 0.0
        total = 0.0
        for z in z_layers:
            tau_z = np.exp(-kappa_e * (height - z) / np.cos(theta))
            total += n_leaf * sigma_leaf * tau_z**2 * dz
        return total

    def _compute_two_way_attenuation(self, veg_medium, polarization) -> float:
        LAI = self._canopy.lai
        kappa_e = 0.5 * LAI
        height = self._canopy.height
        theta = self._geometry.theta_i_rad
        return np.exp(-2.0 * kappa_e * height / np.cos(theta))
