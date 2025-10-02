"""Unit tests for the Mironov soil dielectric model."""

from __future__ import annotations

import math

import numpy as np

from mwrtms.medium import MironovSoilMedium, mironov_permittivity


def test_mironov_permittivity_matches_reference() -> None:
    moisture = 0.25
    clay_fraction = 0.3
    frequency_hz = 5.405e9

    eps = mironov_permittivity(moisture, clay_fraction, frequency_hz)
    assert isinstance(eps, complex)
    # Reference value obtained by evaluating the MATLAB implementation.
    assert math.isclose(eps.real, 11.24896653310607, rel_tol=0.0, abs_tol=1e-9)
    assert math.isclose(eps.imag, 2.5254712378593664, rel_tol=0.0, abs_tol=1e-9)


def test_mironov_medium_uses_helper() -> None:
    medium = MironovSoilMedium(moisture=0.25, clay_fraction=0.3)
    freq = 5.405e9
    eps_helper = mironov_permittivity(0.25, 0.3, freq)
    eps_medium = medium.permittivity(freq)
    assert np.isclose(eps_medium, eps_helper)
    assert medium.physical_temperature() == 290.0
    assert medium.permeability(freq) == complex(1.0)
