from mwrtms.medium import SoilMedium, VegetationMedium


def test_soil_medium_permittivity_cache():
    soil = SoilMedium(
        moisture_m3m3=0.25,
        clay_fraction=0.3,
        sand_fraction=0.5,
        dielectric_model="mironov",
    )
    eps1 = soil.permittivity(5.0e9)
    eps2 = soil.permittivity(5.0e9)
    assert eps1 == eps2


def test_vegetation_medium_permittivity():
    veg = VegetationMedium(gravimetric_moisture=0.6)
    eps = veg.permittivity(5.0e9)
    assert isinstance(eps, complex)
