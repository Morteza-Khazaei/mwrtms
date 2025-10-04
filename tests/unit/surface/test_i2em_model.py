import numpy as np

from mwrtms.core import ElectromagneticWave, PolarizationState, ScatteringGeometry
from mwrtms.medium import HomogeneousMedium
from mwrtms.medium.surface import build_surface_from_statistics
from mwrtms.scattering.surface import I2EMModel


def test_i2em_backscatter_matches_reference():
    wave = ElectromagneticWave(frequency_hz=4.75e9)
    geometry = ScatteringGeometry(theta_i_deg=30.0)
    surface = build_surface_from_statistics(0.009, 0.082)
    model = I2EMModel(wave, geometry, surface, correlation_type="exponential")

    air = HomogeneousMedium(1.0 + 0.0j)
    soil = HomogeneousMedium(15.2 - 2.12j)

    result = model.run(air, soil)

    assert np.isclose(result[PolarizationState.VV], 2.284146950199358e-04, rtol=1e-3)
    assert np.isclose(result[PolarizationState.HH], 1.5226194919292625e-04, rtol=1e-3)
