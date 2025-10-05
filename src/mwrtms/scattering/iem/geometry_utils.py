"""Geometry utilities for IEM family models."""

from __future__ import annotations

import numpy as np

__all__ = [
    "compute_q_vectors",
    "compute_slope_components",
    "compute_spatial_frequency",
]


def compute_q_vectors(
    k: float,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    phi_i: float = 0.0,
) -> tuple[float, float, float, float]:
    """
    Compute incident and scattered wave vector components.
    
    Parameters
    ----------
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
    kx : float
        Incident wave vector x-component
    ky : float
        Incident wave vector y-component
    ksx : float
        Scattered wave vector x-component
    ksy : float
        Scattered wave vector y-component
    
    Notes
    -----
    Convention: z-axis points upward (away from surface)
    - Incident: k_i = k * (sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), -cos(theta_i))
    - Scattered: k_s = k * (sin(theta_s)*cos(phi_s), sin(theta_s)*sin(phi_s), cos(theta_s))
    """
    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    
    kx = k * si * np.cos(phi_i)
    ky = k * si * np.sin(phi_i)
    ksx = k * sis * np.cos(phi_s)
    ksy = k * sis * np.sin(phi_s)
    
    return kx, ky, ksx, ksy


def compute_slope_components(
    theta_i: float,
    theta_s: float,
    phi_s: float,
    cs: float,
    css: float,
    phi_i: float = 0.0,
) -> tuple[float, float]:
    """
    Compute normalized surface slope components.
    
    Parameters
    ----------
    theta_i : float
        Incident angle in radians
    theta_s : float
        Scattered angle in radians
    phi_s : float
        Scattered azimuth angle in radians
    cs : float
        cos(theta_i)
    css : float
        cos(theta_s)
    phi_i : float, optional
        Incident azimuth angle in radians (default: 0.0)
    
    Returns
    -------
    zxx : float
        Normalized slope x-component
    zyy : float
        Normalized slope y-component
    
    Notes
    -----
    From AIEM.m:
    - zxx = -(sin(theta_s)*cos(phi_s) - sin(theta_i)) / (cos(theta_s) + cos(theta_i))
    - zyy = -sin(theta_s)*sin(phi_s) / (cos(theta_s) + cos(theta_i))
    """
    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    csfs = np.cos(phi_s)
    sfs = np.sin(phi_s)
    
    denom = css + cs
    if np.abs(denom) < 1e-10:
        # Grazing angle case
        zxx = 0.0
        zyy = 0.0
    else:
        zxx = -(sis * csfs - si * np.cos(phi_i)) / denom
        zyy = -(sis * sfs - si * np.sin(phi_i)) / denom
    
    return zxx, zyy


def compute_spatial_frequency(
    kl: float,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    phi_i: float = 0.0,
) -> float:
    """
    Compute spatial frequency component K for roughness spectrum.
    
    Parameters
    ----------
    kl : float
        Normalized correlation length (k * L)
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
    K : float
        Spatial frequency magnitude
    
    Notes
    -----
    From AIEM.m:
    K = kl * sqrt((sin(theta_s)*cos(phi_s) - sin(theta_i)*cos(phi_i))^2 
                  + (sin(theta_s)*sin(phi_s) - sin(theta_i)*sin(phi_i))^2)
    """
    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    csfs = np.cos(phi_s)
    sfs = np.sin(phi_s)
    cf = np.cos(phi_i)
    sf = np.sin(phi_i)
    
    K = kl * np.sqrt((sis * csfs - si * cf)**2 + (sis * sfs - si * sf)**2)
    
    return K
