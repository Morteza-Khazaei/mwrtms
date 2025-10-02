from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState
from mwrtms.result import ScatteringResult


def test_facade_soil_backscatter_minimal():
    config = RadarConfigurationFactory.create_monostatic(40.0)
    sigma_vv = mwRTMs.compute_soil_backscatter(
        model="spm",
        radar_config=config,
        frequency_ghz=5.4,
        rms_height_cm=0.5,
        correlation_length_cm=5.0,
        soil_permittivity=10.0 + 1.0j,
        polarization=PolarizationState.VV,
    )
    assert sigma_vv > 0.0


def test_facade_returns_result_for_multiple_channels():
    config = RadarConfigurationFactory.create_monostatic(35.0)
    result = mwRTMs.compute_soil_backscatter(
        model="spm",
        radar_config=config,
        frequency_ghz=3.0,
        rms_height_cm=0.7,
        correlation_length_cm=4.5,
        soil_permittivity=8.0 + 0.9j,
    )
    assert isinstance(result, ScatteringResult)
    assert result[PolarizationState.VV] > 0.0
    assert result[PolarizationState.HH] > 0.0
