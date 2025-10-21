import numpy as np

from mwrtms.scattering.surface.iem.complementary import compute_complementary_vv


def test_complementary_vv_direction_sign_flip() -> None:
    """Ensure downward-propagating branch applies the correct sign to q_i in denominators."""
    u = 0.0
    v = 0.0
    q_mag = 0.8
    Rv = 0.2 + 0.0j
    eps_r = 6.0 + 0.0j
    si = 0.3
    sis = 0.3
    cs = np.sqrt(1.0 - si**2)
    css = np.sqrt(1.0 - sis**2)
    phi = 0.0
    sfs = np.sin(phi)
    csfs = np.cos(phi)

    downward_correct = compute_complementary_vv(
        u,
        v,
        -q_mag,
        -q_mag,
        q_mag,
        Rv,
        eps_r,
        si,
        sis,
        cs,
        css,
        sfs,
        csfs,
        direction=-1,
        is_substrate=False,
    )

    downward_wrong = compute_complementary_vv(
        u,
        v,
        -q_mag,
        -q_mag,
        q_mag,
        Rv,
        eps_r,
        si,
        sis,
        cs,
        css,
        sfs,
        csfs,
        direction=1,  # Missing sign flip (old behaviour)
        is_substrate=False,
    )

    assert np.allclose(downward_correct, -downward_wrong, rtol=1e-12, atol=1e-12)
