import numpy as np

from mwrtms.scattering.surface.iem.multiple_scattering import (
    SurfaceParams,
    _build_propagators,
    _compute_vertical_wavenumbers,
    _evaluate_propagator_at_shifted_point,
    _prepare_geometry_params,
)


def _sample_inputs():
    theta_i = np.deg2rad(35.0)
    theta_s = np.deg2rad(25.0)
    phi_i = 0.0
    phi_s = 0.6
    k = 12.0
    er = 5.0 - 0.8j

    geom = _prepare_geometry_params(theta_i, theta_s, phi_i, phi_s, k)

    U = np.array([[0.1, 0.2], [0.3, 0.4]])
    V = np.array([[0.05, 0.1], [0.15, 0.2]])

    q1, q2 = _compute_vertical_wavenumbers(U, V, k, er)
    propagators = _build_propagators(
        U,
        V,
        q1,
        q2,
        k,
        er,
        geom,
        "vv",
        enable_guardrails=False,
    )

    return geom, k, er, U, V, q1, q2, propagators


def test_propagator_shift_matches_direct_evaluation():
    geom, k, er, U, V, q1, q2, base_props = _sample_inputs()

    def expected_uv(term_index: int):
        norm = np.sqrt(geom.kx**2 + geom.ky**2 + geom.kz**2)
        if norm == 0:
            norm = 1.0
        kx = geom.kx / norm
        ky = geom.ky / norm
        ksx = geom.ksx / norm
        ksy = geom.ksy / norm

        if term_index == 2:
            return (
                -kx - ksx - U,
                -ky - ksy - V,
            )
        if term_index in {3, 4, 5, 12, 13, 14}:
            return (
                -ksx - U,
                -ksy - V,
            )
        if term_index in {6, 7, 8, 9, 10, 11}:
            return (
                -kx - U,
                -ky - V,
            )
        raise AssertionError(f"Unsupported term_index {term_index}")

    for term_index in (2, 3, 6, 9, 12):
        shifted_prop = _evaluate_propagator_at_shifted_point(
            base_props,
            "Fp",
            term_index,
            U,
            V,
            k,
            er,
            geom,
            "vv",
            enable_guardrails=False,
        )

        U_expected, V_expected = expected_uv(term_index)
        q1_expected, q2_expected = _compute_vertical_wavenumbers(
            U_expected,
            V_expected,
            k,
            er,
        )
        expected_props = _build_propagators(
            U_expected,
            V_expected,
            q1_expected,
            q2_expected,
            k,
            er,
            geom,
            "vv",
            enable_guardrails=False,
        )

        np.testing.assert_allclose(
            shifted_prop,
            expected_props["Fp"],
            rtol=1e-9,
            atol=1e-12,
        )
