import numpy as np

from mwrtms.core import ElectromagneticWave, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.core import RadarConfigurationFactory
from mwrtms.scattering.surface import SPMModel
from mwrtms.core.polarization import PolarizationState
from mwrtms.factory import ScatteringModelFactory


def test_spm_cross_pol_returns_positive_power() -> None:
    surface = build_surface_from_statistics(
        rms_height_m=0.006,
        correlation_length_m=0.03,
        correlation_type="exponential",
    )

    wave = ElectromagneticWave(frequency_hz=5.4e9)
    geometry = ScatteringGeometry(theta_i_deg=40.0)

    air = HomogeneousMedium(1.0 + 0j)
    soil = HomogeneousMedium(5.0 + 1.0j)

    model = SPMModel(wave, geometry, surface)

    result = model.run(air, soil)
    assert result[PolarizationState.HV] > 0.0
    assert result[PolarizationState.HH] > 0.0
    assert result[PolarizationState.VV] > result[PolarizationState.HV]



def test_compute_with_config_validates_geometry() -> None:
    surface = build_surface_from_statistics(0.004, 0.02)
    wave = ElectromagneticWave(frequency_hz=2.0e9)
    air = HomogeneousMedium(1.0 + 0.0j)
    soil = HomogeneousMedium(7.0 + 1.2j)

    config_ref = RadarConfigurationFactory.create_monostatic(30.0)
    config_other = RadarConfigurationFactory.create_monostatic(35.0)

    model = ScatteringModelFactory.create_with_radar_config(
        "spm",
        config=config_ref,
        wave=wave,
        surface=surface,
    )

    # Matching configuration succeeds
    result = model.compute_with_config(air, soil, None, config_ref)
    assert result[PolarizationState.VV] > 0.0

    # Mismatched geometry raises
    import pytest

    with pytest.raises(ValueError):
        model.compute_with_config(air, soil, PolarizationState.VV, config_other)
