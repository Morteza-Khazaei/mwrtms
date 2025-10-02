from mwrtms import mwRTMs


def test_facade_soil_backscatter_minimal():
    result = mwRTMs.compute_soil_backscatter(
        model="spm",
        frequency_ghz=5.4,
        incident_angle_deg=40.0,
        rms_height_cm=0.5,
        correlation_length_cm=5.0,
        soil_moisture=0.2,
    )
    for pol in ("hh", "vv", "hv"):
        assert pol in result
