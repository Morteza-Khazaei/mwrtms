from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState


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
