import numpy as np

from mwrtms.core import ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import Medium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface import SPMModel
from mwrtms.core.polarization import PolarizationState


class _ConstantMedium(Medium):
    def __init__(self, permittivity: complex) -> None:
        super().__init__(temperature_k=293.15)
        self._eps = complex(permittivity)

    def permittivity(self, frequency_hz: float) -> complex:
        return self._eps


def test_spm_cross_pol_returns_positive_power() -> None:
    surface = build_surface_from_statistics(
        rms_height_m=0.006,
        correlation_length_m=0.03,
        correlation_type="exponential",
    )

    wave = ElectromagneticWave(frequency_hz=5.4e9)
    geometry = ScatteringGeometry(theta_i_deg=40.0)

    air = _ConstantMedium(1.0 + 0j)
    soil = _ConstantMedium(5.0 + 1.0j)

    model = SPMModel(wave, geometry, surface)

    sigma_hv = model.compute(air, soil, PolarizationState.HV)
    assert sigma_hv > 0.0

    sigma_hh = model.compute(air, soil, PolarizationState.HH)
    sigma_vv = model.compute(air, soil, PolarizationState.VV)

    assert sigma_hh > 0.0 and sigma_vv > 0.0
    assert sigma_hv < sigma_vv
