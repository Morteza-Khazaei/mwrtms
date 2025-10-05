import os
import pytest
import numpy as np


@pytest.fixture(autouse=True, scope="session")
def plan_a_fast_profile():
    """
    Enforce the fast Numba profile (Plan A) for I2EM during test session.

    Settings chosen to replicate the earlier fast and accurate behavior:
      - Use Numba backend
      - Gauss–Legendre quadrature (fewer nodes needed for accuracy)
      - Fast kernel (fastmath=True)
      - Modest node counts (r=81, phi=97)
      - No adaptive refinement by default

    Domain alignment for the fast profile is handled in code; ensure tests
    do not override with heavier ranges unless explicitly required.
    """
    os.environ["MWRTMS_I2EM_USE_NUMBA"] = "1"
    os.environ["MWRTMS_I2EM_QUADRATURE"] = "gauss"
    os.environ["MWRTMS_I2EM_NUMBA_FASTMATH"] = "1"
    os.environ["MWRTMS_I2EM_R_POINTS"] = "81"
    os.environ["MWRTMS_I2EM_PHI_POINTS"] = "97"
    # Ensure adaptive refinement is disabled for speed in regression tests
    os.environ.pop("MWRTMS_I2EM_ADAPTIVE", None)
    # Optional: tolerance if enabled elsewhere; ensure default is not set
    os.environ.pop("MWRTMS_I2EM_ADAPT_TOL", None)

    # Ensure the correlation model used in tests aligns with NMM3D (exponential).
    # If your tests explicitly set correlation, this line is harmless.
    os.environ.setdefault("MWRTMS_I2EM_CORRELATION", "exponential")

    yield


@pytest.fixture(autouse=True, scope="session")
def deterministic_seed():
    """Seed numpy RNG for reproducible test sampling, if any."""
    np.random.seed(12345)
    yield
