"""Kirchhoff field coefficients for AIEM."""

from __future__ import annotations

import numpy as np

__all__ = ["compute_kirchhoff_coefficients"]


def compute_kirchhoff_coefficients(
    Rv: complex,
    Rh: complex,
    k: float,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    phi_i: float = 0.0,
) -> tuple[complex, complex, complex, complex]:
    """
    Compute Kirchhoff field coefficients for all polarizations.
    
    Parameters
    ----------
    Rv : complex
        Vertical polarization reflection coefficient (transition-adjusted)
    Rh : complex
        Horizontal polarization reflection coefficient (transition-adjusted)
    k : float
        Wavenumber (2π/λ)
    theta_i : float
        Incident angle in radians
    theta_s : float
        Scattered angle in radians
    phi_s : float
        Scattered azimuth angle in radians
    phi_i : float, optional
        Incident azimuth angle in radians (default: 0.0)
    
    Returns
    -------
    fvv : complex
        VV polarization Kirchhoff field coefficient
    fhh : complex
        HH polarization Kirchhoff field coefficient
    fhv : complex
        HV polarization Kirchhoff field coefficient
    fvh : complex
        VH polarization Kirchhoff field coefficient
    
    Notes
    -----
    From AIEM.m, the Kirchhoff field coefficients are computed using
    the local surface slope components and reflection coefficients.
    
    The formulation accounts for:
    - Surface slope effects (zxx, zyy)
    - Polarization coupling
    - Geometric factors from incident and scattered directions
    """
    # Trigonometric quantities
    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    cs = np.cos(theta_i)
    css = np.cos(theta_s)
    cf = np.cos(phi_i)
    csfs = np.cos(phi_s)
    sfs = np.sin(phi_s)
    sf = np.sin(phi_i)
    
    cs2 = cs * cs
    si2 = si * si
    
    # Slope components
    denom = css + cs
    if np.abs(denom) < 1e-10:
        # Grazing angle - return zeros
        return 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j
    
    zxx = -(sis * csfs - si * cf) / denom
    zyy = -(sis * sfs - si * sf) / denom
    
    # Distance factor
    d2 = np.sqrt((zxx * cs - si * cf) ** 2 + (zyy * cs) ** 2 + (cs + si * zxx) ** 2)
    if np.abs(d2) < 1e-10:
        d2 = 1e-10
    
    # Field components for HH polarization
    hsnv = -(cs * csfs + si * (zxx * csfs + zyy * sfs))
    vsnh = css * csfs - zxx * sis
    hsnh = -sfs
    vsnv = zyy * cs * sis + css * (zyy * csfs * si - (cs + zxx * si) * sfs)
    
    # Tangential and normal components
    hsnt = (-(cs2 + si2) * sfs * (-si * cf + cs * zxx) + 
            csfs * (cs + si * zxx) * zyy + si * sfs * (zyy ** 2)) / d2
    
    hsnd = (-(cs + si * zxx) * (-csfs * si * cf + cs * csfs * zxx + cs * sfs * zyy)) / d2
    
    vsnt = ((cs2 + si2) * (-si * cf + cs * zxx) * (csfs * css - sis * zxx) + 
            css * sfs * (cs + si * zxx) * zyy - 
            (csfs * css * si * cf + cs * sis) * (zyy ** 2)) / d2
    
    vsnd = (-(cs + si * zxx) * (si * sis * zyy - 
            css * (si * sfs - cs * sfs * zxx + cs * csfs * zyy))) / d2
    
    # Kirchhoff field coefficients
    fhh = ((1.0 - Rh) * hsnv + (1.0 + Rh) * vsnh - 
           (hsnt + vsnd) * (Rh + Rv) * (zyy / d2))
    
    fvv = -((1.0 - Rv) * hsnv + (1.0 + Rv) * vsnh) + \
          (hsnt + vsnd) * (Rh + Rv) * (zyy / d2)
    
    fhv = (-(1.0 + Rv) * hsnh + (1.0 - Rv) * vsnv + 
           (hsnd - vsnt) * (Rh + Rv) * (zyy / d2))
    
    fvh = (-(1.0 + Rh) * hsnh + (1.0 - Rh) * vsnv + 
           (hsnd - vsnt) * (Rh + Rv) * (zyy / d2))
    
    return fvv, fhh, fhv, fvh
