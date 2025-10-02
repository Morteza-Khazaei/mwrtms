"""Tests demonstrating the four OOP pillars implemented in mwRTMs."""

from __future__ import annotations

import pathlib
import sys

SYS_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SYS_ROOT / "src"))

import math

import pytest

from mwrtms import (
    AIEMModel,
    ElectromagneticWave,
    ScatteringGeometry,
    SurfaceScattering,
    build_surface_from_statistics,
    RadarConfigurationFactory,
    PolarizationState,
)
from mwrtms.factory import ScatteringModelFactory
from mwrtms.facade import mwRTMs
from mwrtms.result import ScatteringResult
from mwrtms.medium import HomogeneousMedium


@pytest.fixture()
def simple_setup():
    wave = ElectromagneticWave(5.4e9)
    geometry = ScatteringGeometry(theta_i_deg=40.0)
    surface = build_surface_from_statistics(0.01, 0.05, correlation_type="exponential")
    air = HomogeneousMedium(1.0 + 0.0j)
    soil = HomogeneousMedium(15 - 2j)
    return wave, geometry, surface, air, soil


def test_encapsulation(simple_setup):
    wave, geometry, surface, *_ = simple_setup
    with pytest.raises(AttributeError):
        wave._frequency_hz = 10e9  # type: ignore[attr-defined]
    assert math.isclose(geometry.theta_i, math.radians(40.0))
    assert pytest.approx(surface.rms_height()) == 0.01


def test_inheritance_polymorphism(simple_setup):
    wave, geometry, surface, air, soil = simple_setup
    model = AIEMModel(wave, geometry, surface)
    assert isinstance(model, SurfaceScattering)
    result = model.run(air, soil)
    assert isinstance(result, ScatteringResult)
    assert result[PolarizationState.VV] >= 0.0

    factory_model = ScatteringModelFactory.create(
        "aiem", wave=wave, geometry=geometry, surface=surface
    )
    result_factory = factory_model.run(air, soil)
    assert result_factory[PolarizationState.VV] >= 0.0


def test_abstraction_facade(simple_setup):
    config = RadarConfigurationFactory.create_monostatic(40.0)
    sigma_vv = mwRTMs.compute_soil_backscatter(
        model="aiem",
        radar_config=config,
        frequency_ghz=5.4,
        rms_height_cm=1.0,
        correlation_length_cm=5.0,
        soil_permittivity=12.0 + 2.0j,
        polarization=PolarizationState.VV,
    )
    assert sigma_vv > 0.0
