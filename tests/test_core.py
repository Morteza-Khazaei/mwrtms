import math

from mwrtms.core import ElectromagneticWave, ScatteringGeometry


def test_wave_properties():
    wave = ElectromagneticWave(10.0e9)
    assert math.isclose(wave.frequency_ghz, 10.0)
    assert wave.wavelength > 0
    assert wave.wavenumber > 0


def test_geometry_backscatter():
    geometry = ScatteringGeometry(theta_i_deg=40.0, theta_s_deg=40.0, phi_s_deg=180.0)
    assert geometry.is_backscatter
    vectors = geometry.wave_vectors(geometry.theta_i_rad + 1.0)  # arbitrary k
    assert set(vectors.keys()) == {"ki", "ks"}
