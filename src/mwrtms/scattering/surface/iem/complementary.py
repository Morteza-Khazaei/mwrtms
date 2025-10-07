"""Complementary field coefficients for AIEM multiple scattering term."""

from __future__ import annotations

import numpy as np

__all__ = [
    "compute_expal",
    "compute_complementary_vv",
    "compute_complementary_hh",
    "compute_complementary_hv",
    "compute_complementary_vh",
]


def compute_expal(q: complex, ks: float, cs: float, css: float) -> complex:
    """
    Compute exponential factor for complementary field coefficients.
    
    Parameters
    ----------
    q : complex
        Wave vector component
    ks : float
        Normalized RMS height (k * sigma)
    cs : float
        cos(theta_i)
    css : float
        cos(theta_s)
    
    Returns
    -------
    expal : complex
        Exponential factor exp(-ks^2 * (q^2 - q*(css - cs)))
    
    Notes
    -----
    From AIEM.m:
    expal = exp(-ks^2 * (q^2 - q*(css - cs)))
    """
    ks2 = ks * ks
    return np.exp(-ks2 * (q * q - q * (css - cs)))


def _compute_slope_and_field_components(
    u: float,
    v: float,
    q: complex,
    qslp: complex,
    si: float,
    sis: float,
    cs: float,
    css: float,
    sfs: float,
    csfs: float,
) -> tuple:
    """
    Compute slope components and intermediate field quantities.
    
    Returns
    -------
    zx, zy, zxp, zyp : complex
        Slope components
    """
    kxu = si + u
    ksxu = sis * csfs + u
    kyv = v
    ksyv = sis * sfs + v
    
    # Slope components at scattered angle
    # CORRECTED: Use complex magnitude, not real part magnitude
    if np.abs(css - qslp) < 1e-10:
        zx = 0.0 + 0.0j
        zy = 0.0 + 0.0j
    else:
        zx = -ksxu / (css - qslp)
        zy = -ksyv / (css - qslp)
    
    # Slope components at incident angle
    # CORRECTED: Use complex magnitude, not real part magnitude
    if np.abs(cs + qslp) < 1e-10:
        zxp = 0.0 + 0.0j
        zyp = 0.0 + 0.0j
    else:
        zxp = kxu / (cs + qslp)
        zyp = kyv / (cs + qslp)
    
    return zx, zy, zxp, zyp


def compute_complementary_vv(
    u: float,
    v: float,
    q: complex,
    qslp: complex,
    qfix: complex,
    Rv: complex,
    eps_r: complex,
    si: float,
    sis: float,
    cs: float,
    css: float,
    sfs: float,
    csfs: float,
    is_substrate: bool = False,
) -> complex:
    """
    Compute complementary field coefficient for VV polarization.
    
    Parameters
    ----------
    u, v : float
        Integration variables
    q : complex
        Wave vector component
    qslp : complex
        Slope-dependent wave vector
    qfix : complex
        Fixed wave vector component
    Rv : complex
        Vertical polarization reflection coefficient
    eps_r : complex
        Relative dielectric constant (used only for substrate-side)
    si, sis : float
        sin(theta_i), sin(theta_s)
    cs, css : float
        cos(theta_i), cos(theta_s)
    sfs, csfs : float
        sin(phi_s), cos(phi_s)
    is_substrate : bool
        If True, compute substrate-side (Fbvv), else air-side (Favv)
    
    Returns
    -------
    Fvv : complex
        Complementary field coefficient for VV polarization
    
    Notes
    -----
    Implements favv() and fbvv() from AIEM.m
    """
    zx, zy, zxp, zyp = _compute_slope_and_field_components(
        u, v, q, qslp, si, sis, cs, css, sfs, csfs
    )
    
    # Field components (same for air and substrate)
    c1 = -csfs * (-1.0 - zx * zxp) + sfs * zxp * zy
    
    c2 = (-csfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - 
                   cs * v * zyp - si * v * zx * zyp) +
          sfs * (cs * u * zy + si * u * zxp * zy + q * si * zyp - 
                 cs * u * zyp + si * v * zy * zyp))
    
    c3 = (-csfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) +
          sfs * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy))
    
    c4 = (-css * sfs * (-si * zyp + cs * zx * zyp) - 
          csfs * css * (-cs - si * zxp - cs * zy * zyp) +
          sis * (-cs * zx - si * zx * zxp - si * zy * zyp))
    
    c5 = (-css * sfs * (-v * zx + v * zxp) - 
          csfs * css * (q + u * zxp + v * zy) +
          sis * (q * zx + u * zx * zxp + v * zxp * zy))
    
    c6 = (-css * sfs * (-u * zyp + q * zx * zyp) - 
          csfs * css * (v * zyp - q * zy * zyp) +
          sis * (v * zx * zyp - u * zy * zyp))
    
    # Reflection coefficient factors
    rpv = 1.0 + Rv
    rmv = 1.0 - Rv
    av = rpv / qfix
    bv = rmv / qfix
    
    if not is_substrate:
        # Air-side: favv
        Fvv = bv * (-rpv * c1 + rmv * c2 + rpv * c3) + av * (rmv * c4 + rpv * c5 + rmv * c6)
    else:
        # Substrate-side: fbvv
        Fvv = (av * (rpv * c1 - rmv * c2 - rpv * c3 / eps_r) - 
               bv * (rmv * c4 * eps_r + rpv * c5 + rmv * c6))
    
    return Fvv


def compute_complementary_hh(
    u: float,
    v: float,
    q: complex,
    qslp: complex,
    qfix: complex,
    Rh: complex,
    eps_r: complex,
    si: float,
    sis: float,
    cs: float,
    css: float,
    sfs: float,
    csfs: float,
    is_substrate: bool = False,
) -> complex:
    """
    Compute complementary field coefficient for HH polarization.
    
    Parameters are the same as compute_complementary_vv, but using Rh instead of Rv.
    
    Notes
    -----
    Implements fahh() and fbhh() from AIEM.m
    """
    zx, zy, zxp, zyp = _compute_slope_and_field_components(
        u, v, q, qslp, si, sis, cs, css, sfs, csfs
    )
    
    # Field components (same as VV)
    c1 = -csfs * (-1.0 - zx * zxp) + sfs * zxp * zy
    
    c2 = (-csfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - 
                   cs * v * zyp - si * v * zx * zyp) +
          sfs * (cs * u * zy + si * u * zxp * zy + q * si * zyp - 
                 cs * u * zyp + si * v * zy * zyp))
    
    c3 = (-csfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) +
          sfs * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy))
    
    c4 = (-css * sfs * (-si * zyp + cs * zx * zyp) - 
          csfs * css * (-cs - si * zxp - cs * zy * zyp) +
          sis * (-cs * zx - si * zx * zxp - si * zy * zyp))
    
    c5 = (-css * sfs * (-v * zx + v * zxp) - 
          csfs * css * (q + u * zxp + v * zy) +
          sis * (q * zx + u * zx * zxp + v * zxp * zy))
    
    c6 = (-css * sfs * (-u * zyp + q * zx * zyp) - 
          csfs * css * (v * zyp - q * zy * zyp) +
          sis * (v * zx * zyp - u * zy * zyp))
    
    # Reflection coefficient factors
    rph = 1.0 + Rh
    rmh = 1.0 - Rh
    ah = rph / qfix
    bh = rmh / qfix
    
    if not is_substrate:
        # Air-side: fahh
        Fhh = -bh * (-rph * c1 + rmh * c2 + rph * c3) - ah * (rmh * c4 + rph * c5 + rmh * c6)
    else:
        # Substrate-side: fbhh
        Fhh = (ah * (-rph * c1 * eps_r + rmh * c2 + rph * c3) + 
               bh * (rmh * c4 + rph * c5 + rmh * c6 / eps_r))
    
    return Fhh


def compute_complementary_hv(
    u: float,
    v: float,
    q: complex,
    qslp: complex,
    qfix: complex,
    Rhv: complex,
    eps_r: complex,
    si: float,
    sis: float,
    cs: float,
    css: float,
    sfs: float,
    csfs: float,
    is_substrate: bool = False,
) -> complex:
    """
    Compute complementary field coefficient for HV polarization.
    
    Notes
    -----
    Implements fahv() and fbhv() from AIEM.m
    """
    zx, zy, zxp, zyp = _compute_slope_and_field_components(
        u, v, q, qslp, si, sis, cs, css, sfs, csfs
    )
    
    # Field components for cross-polarization
    b1 = (-css * sfs * (-1.0 - zx * zxp) - sis * zy - csfs * css * zxp * zy)
    
    b2 = (-css * sfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - 
                        cs * v * zyp - si * v * zx * zyp) +
          sis * (-cs * q * zy - q * si * zxp * zy + q * si * zx * zyp - 
                 cs * u * zx * zyp - cs * v * zy * zyp) -
          csfs * css * (cs * u * zy + si * u * zxp * zy + q * si * zyp - 
                        cs * u * zyp + si * v * zy * zyp))
    
    b3 = (-css * sfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) -
          csfs * css * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy) +
          sis * (-si * v * zx + cs * v * zx * zxp + si * u * zy - cs * u * zxp * zy))
    
    b4 = -csfs * (-si * zyp + cs * zx * zyp) + sfs * (-cs - si * zxp - cs * zy * zyp)
    
    b5 = -csfs * (-v * zx + v * zxp) + sfs * (q + u * zxp + v * zy)
    
    b6 = -csfs * (-u * zyp + q * zx * zyp) + sfs * (v * zyp - q * zy * zyp)
    
    # Reflection coefficient factors
    rp = 1.0 + Rhv
    rm = 1.0 - Rhv
    a = rp / qfix
    b = rm / qfix
    
    if not is_substrate:
        # Air-side: fahv
        Fhv = b * (rp * b1 - rm * b2 - rp * b3) + a * (rm * b4 + rp * b5 + rm * b6)
    else:
        # Substrate-side: fbhv
        Fhv = (a * (-rp * b1 + rm * b2 + rp * b3 / eps_r) - 
               b * (rm * b4 * eps_r + rp * b5 + rm * b6))
    
    return Fhv


def compute_complementary_vh(
    u: float,
    v: float,
    q: complex,
    qslp: complex,
    qfix: complex,
    Rhv: complex,
    eps_r: complex,
    si: float,
    sis: float,
    cs: float,
    css: float,
    sfs: float,
    csfs: float,
    is_substrate: bool = False,
) -> complex:
    """
    Compute complementary field coefficient for VH polarization.
    
    Notes
    -----
    Implements favh() and fbvh() from AIEM.m
    """
    zx, zy, zxp, zyp = _compute_slope_and_field_components(
        u, v, q, qslp, si, sis, cs, css, sfs, csfs
    )
    
    # Field components (same as HV)
    b1 = (-css * sfs * (-1.0 - zx * zxp) - sis * zy - csfs * css * zxp * zy)
    
    b2 = (-css * sfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - 
                        cs * v * zyp - si * v * zx * zyp) +
          sis * (-cs * q * zy - q * si * zxp * zy + q * si * zx * zyp - 
                 cs * u * zx * zyp - cs * v * zy * zyp) -
          csfs * css * (cs * u * zy + si * u * zxp * zy + q * si * zyp - 
                        cs * u * zyp + si * v * zy * zyp))
    
    b3 = (-css * sfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) -
          csfs * css * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy) +
          sis * (-si * v * zx + cs * v * zx * zxp + si * u * zy - cs * u * zxp * zy))
    
    b4 = -csfs * (-si * zyp + cs * zx * zyp) + sfs * (-cs - si * zxp - cs * zy * zyp)
    
    b5 = -csfs * (-v * zx + v * zxp) + sfs * (q + u * zxp + v * zy)
    
    b6 = -csfs * (-u * zyp + q * zx * zyp) + sfs * (v * zyp - q * zy * zyp)
    
    # Reflection coefficient factors
    rp = 1.0 + Rhv
    rm = 1.0 - Rhv
    a = rp / qfix
    b = rm / qfix
    
    if not is_substrate:
        # Air-side: favh
        Fvh = b * (rp * b4 + rm * b5 + rp * b6) - a * (-rm * b1 + rp * b2 + rm * b3)
    else:
        # Substrate-side: fbvh
        Fvh = (-a * (rp * b4 + rm * b5 + rp * b6 / eps_r) + 
               b * (-rm * b1 * eps_r + rp * b2 + rm * b3))
    
    return Fvh
