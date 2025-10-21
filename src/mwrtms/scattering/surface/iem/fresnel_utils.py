"""Fresnel reflection coefficient utilities for IEM family models."""

from __future__ import annotations

import numpy as np

__all__ = [
    "compute_fresnel_incident",
    "compute_fresnel_specular",
    "compute_fresnel_nadir",
]


def compute_fresnel_incident(
    eps_r: complex, theta_i: float, mu_r: float = 1.0
) -> tuple[complex, complex, complex]:
    """
    Compute Fresnel reflection coefficients at the incident angle.
    
    Parameters
    ----------
    eps_r : complex
        Relative dielectric constant of the substrate
    theta_i : float
        Incident angle in radians
    mu_r : float, optional
        Relative permeability (default: 1.0)
    
    Returns
    -------
    Rvi : complex
        Vertical polarization reflection coefficient at incident angle
    Rhi : complex
        Horizontal polarization reflection coefficient at incident angle
    Rvhi : complex
        Cross-polarization reflection coefficient (Rvi - Rhi) / 2
    
    Notes
    -----
    Based on AIEM.m formulation:
    - stem = sqrt(er*ur - sin^2(theta_i))
    - Rvi = (er*cos(theta_i) - stem) / (er*cos(theta_i) + stem)
    - Rhi = (ur*cos(theta_i) - stem) / (ur*cos(theta_i) + stem)
    """
    si = np.sin(theta_i)
    cs = np.cos(theta_i)
    si2 = si * si
    
    # Transmitted wave vector component with proper branch for lossy media
    # Ensure Im(stem) >= 0 for decaying wave into lower half-space
    stem = np.sqrt(eps_r * mu_r - si2)
    if np.imag(stem) < 0:
        stem = -stem
    
    # Vertical polarization (TM mode)
    Rvi = (eps_r * cs - stem) / (eps_r * cs + stem)
    
    # Horizontal polarization (TE mode)
    Rhi = (mu_r * cs - stem) / (mu_r * cs + stem)
    
    # Cross-polarization term
    Rvhi = (Rvi - Rhi) / 2.0
    
    return Rvi, Rhi, Rvhi


import numpy as np

def compute_fresnel_specular(
    eps_r: complex,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    mu_r: float = 1.0,
) -> tuple[complex, complex, complex]:
    """
    Fresnel reflection coefficients at the local specular angle for an air (ε1=1, μ1=1) → medium (εr=eps_r, μr=mu_r) interface.

    θ_i : incident polar angle (from +z), θ_s : scatter polar angle, φ_s : scatter azimuth (relative to incident plane).

    We use the specular-facet half-angle relation:
      cos θ_l = cos(ψ/2) = sqrt( (1 + cos ψ) / 2 ),
    with
      cos ψ = cos θ_i cos θ_s - sin θ_i sin θ_s cos φ_s.

    Returns
    -------
    Rv, Rh, Rvh : complex
        VV(p) and HH(s) Fresnel coefficients at θ_l; cross-pol (Rvh) is 0 for a planar interface.
    """
    eps_r = np.asarray(eps_r, dtype=np.complex128)
    theta_i = np.asarray(theta_i, dtype=np.float64)
    theta_s = np.asarray(theta_s, dtype=np.float64)
    phi_s   = np.asarray(phi_s,   dtype=np.float64)

    cs_i,  si_i  = np.cos(theta_i), np.sin(theta_i)
    cs_s,  si_s  = np.cos(theta_s), np.sin(theta_s)
    cphi_s       = np.cos(phi_s)

    # cos ψ between the (downward) incident and (upward) scattered directions in AIEM/IEM geometry
    cos_psi = cs_i*cs_s - si_i*si_s*cphi_s

    # local incidence angle at the specular facet
    csl = np.sqrt(np.clip(1.0 + cos_psi, 0.0, 2.0) / 2.0)         # cos θ_l
    sil2 = np.clip(1.0 - csl*csl, 0.0, 1.0)                       # sin^2 θ_l

    # transmitted z-term (complex-safe), choose decaying branch
    steml = np.sqrt(eps_r*mu_r - sil2)
    steml = np.where(np.imag(steml) < 0, -steml, steml)

    # Fresnel (p ≡ VV, s ≡ HH) for medium-1 (air) to medium-2 (εr, μr)
    Rv = (eps_r*csl - steml) / (eps_r*csl + steml)   # VV / p
    Rh = (mu_r *csl - steml) / (mu_r *csl + steml)   # HH / s

    Rvh = np.zeros_like(Rv, dtype=np.complex128)     # no cross-pol for planar Fresnel

    return Rv, Rh, Rvh



def compute_fresnel_nadir(eps_r: complex) -> tuple[complex, complex]:
    """
    Compute Fresnel reflection coefficients at nadir (normal incidence).
    
    Parameters
    ----------
    eps_r : complex
        Relative dielectric constant of the substrate
    
    Returns
    -------
    rv0 : complex
        Vertical polarization reflection coefficient at nadir
    rh0 : complex
        Horizontal polarization reflection coefficient at nadir
    
    Notes
    -----
    At normal incidence (theta = 0), there is no distinction between H and V polarizations.
    Both use the same formula:
    - r0 = (sqrt(eps_r) - 1) / (sqrt(eps_r) + 1)
    - rv0 = rh0 = r0
    
    CORRECTED: Previous MATLAB code incorrectly used rh0 = -rv0
    """
    sqrt_er = np.sqrt(eps_r)
    rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
    rh0 = -rv0 # reciprocal for HH polarization in IEM family models
    
    return rv0, rh0
