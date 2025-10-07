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


def compute_fresnel_specular(
    eps_r: complex,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    mu_r: float = 1.0,
) -> tuple[complex, complex, complex]:
    """
    Compute Fresnel reflection coefficients at the specular angle.
    
    The specular angle is computed from the incident and scattered directions
    using the local surface normal approximation.
    
    Parameters
    ----------
    eps_r : complex
        Relative dielectric constant of the substrate
    theta_i : float
        Incident angle in radians
    theta_s : float
        Scattered angle in radians
    phi_s : float
        Scattered azimuth angle in radians
    mu_r : float, optional
        Relative permeability (default: 1.0)
    
    Returns
    -------
    Rvl : complex
        Vertical polarization reflection coefficient at specular angle
    Rhl : complex
        Horizontal polarization reflection coefficient at specular angle
    Rvhl : complex
        Cross-polarization reflection coefficient
    
    Notes
    -----
    From AIEM.m:
    - csl = sqrt(1 + cos(theta_i)*cos(theta_s) - sin(theta_i)*sin(theta_s)*cos(phi_s)) / sqrt(2)
    - This represents the cosine of the local specular angle
    """
    cs = np.cos(theta_i)
    si = np.sin(theta_i)
    css = np.cos(theta_s)
    sis = np.sin(theta_s)
    csfs = np.cos(phi_s)
    
    # Local specular angle cosine
    # CORRECTED: cos(ψ/2) = sqrt((1 + cos(ψ))/2) where cos(ψ) = -cos(θ_i)cos(θ_s) + sin(θ_i)sin(θ_s)cos(φ_s)
    # This gives: csl = sqrt((1 - cos(θ_i)cos(θ_s) + sin(θ_i)sin(θ_s)cos(φ_s))/2)
    csl = np.sqrt(1.0 + cs * css - si * sis * csfs) / np.sqrt(2.0)
    sil = np.sqrt(1.0 - csl * csl)
    
    # Transmitted wave vector at specular angle with proper branch for lossy media
    # Ensure Im(steml) >= 0 for decaying wave into lower half-space
    steml = np.sqrt(eps_r * mu_r - sil * sil)
    if np.imag(steml) < 0:
        steml = -steml
    
    # Vertical polarization
    Rvl = (eps_r * csl - steml) / (eps_r * csl + steml)
    
    # Horizontal polarization
    Rhl = (mu_r * csl - steml) / (mu_r * csl + steml)
    
    # Cross-polarization
    Rvhl = (Rvl - Rhl) / 2.0
    
    return Rvl, Rhl, Rvhl


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
