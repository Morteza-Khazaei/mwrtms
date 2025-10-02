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
    ExponentialCorrelation,
    IsotropicMedium,
    ScatteringGeometry,
    SurfaceRoughness,
)
from mwrtms.core.polarization import PolarizationState
from mwrtms.factory import ScatteringModelFactory
from mwrtms.facade import mwRTMsFacade
from mwrtms.scattering.surface.base import SurfaceScattering


@pytest.fixture()
def simple_setup():
    wave = ElectromagneticWave(5.4e9)
    geometry = ScatteringGeometry(theta_i_deg=40.0)
    roughness = SurfaceRoughness(0.01, 0.05, ExponentialCorrelation(0.05))
    air = IsotropicMedium(1.0, 290.0)
    soil = IsotropicMedium(15 - 2j, 290.0)
    return wave, geometry, roughness, air, soil


def test_encapsulation(simple_setup):
    wave, geometry, roughness, air, soil = simple_setup
    with pytest.raises(AttributeError):
        wave._frequency_hz = 10e9  # type: ignore[attr-defined]
    assert math.isclose(geometry.theta_i, math.radians(40.0))
    assert pytest.approx(roughness.rms_height) == 0.01


def test_inheritance_polymorphism(simple_setup):
    wave, geometry, roughness, air, soil = simple_setup
    model = AIEMModel(wave, geometry, roughness)
    assert isinstance(model, SurfaceScattering)
    sigma = model.compute(air, soil, PolarizationState.VV)
    assert sigma >= 0.0

    factory_model = ScatteringModelFactory.create(
        "i2em", wave=wave, geometry=geometry, surface_roughness=roughness
    )
    sigma_factory = factory_model.compute(air, soil, PolarizationState.VV)
    assert sigma_factory >= 0.0


def test_abstraction_facade(simple_setup):
    _, _, _, _, _ = simple_setup
    sigma = mwRTMsFacade.compute_soil_backscatter(
        model="iem",
        frequency_ghz=5.4,
        incident_angle_deg=40.0,
        rms_height_cm=1.0,
        correlation_length_cm=5.0,
        soil_permittivity=15 - 2j,
    )
    assert sigma >= 0.0
