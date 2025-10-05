"""Transition function for AIEM reflection coefficients."""

from __future__ import annotations

import math
import numpy as np

__all__ = ["compute_transition_function"]


def compute_transition_function(
    eps_r: complex,
    theta_i: float,
    ks: float,
    cs: float,
    spectral_weights: np.ndarray,
    n_terms: int,
) -> tuple[float, float]:
    """
    Compute AIEM transition function for V and H polarizations.
    
    The transition function smoothly interpolates between the Fresnel
    reflection coefficients at the incident angle and the specular angle,
    accounting for surface roughness effects.
    
    Parameters
    ----------
    eps_r : complex
        Relative dielectric constant of the substrate
    theta_i : float
        Incident angle in radians
    ks : float
        Normalized RMS height (k * sigma)
    cs : float
        cos(theta_i)
    spectral_weights : np.ndarray
        Surface roughness spectrum weights W_n
    n_terms : int
        Number of spectral terms
    
    Returns
    -------
    Tfv : float
        Transition function for vertical polarization
    Tfh : float
        Transition function for horizontal polarization
    
    Notes
    -----
    From AIEM.m (T.D. Wu & A.K. Fung formulation):
    
    1. Compute nadir reflection coefficients rv0, rh0
    2. Compute transition factors Ft_v, Ft_h
    3. Compute shadowing terms St, St0
    4. Transition function: Tf = 1 - St/St0
    
    The transition function ensures smooth behavior across the full
    range of incidence angles and surface roughness conditions.
    
    References
    ----------
    Wu, T. D., & Fung, A. K. (1992). A transition model for the reflection
    coefficient in surface scattering. IEEE TGRS, 30(4), 856-860.
    """
    si = np.sin(theta_i)
    s2 = si * si
    
    # Nadir reflection coefficients
    sqrt_er = np.sqrt(eps_r)
    rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
    rh0 = -rv0
    
    # Transmitted wave vector component
    rt = np.sqrt(eps_r - s2)
    
    # Transition factors
    denom_v = cs * rt
    if np.abs(denom_v) < 1e-10:
        Ftv = 0.0
    else:
        Ftv = 8.0 * (rv0 ** 2) * s2 * (cs + rt) / denom_v
    
    denom_h = cs * rt
    if np.abs(denom_h) < 1e-10:
        Fth = 0.0
    else:
        Fth = -8.0 * (rh0 ** 2) * s2 * (cs + rt) / denom_h
    
    # Shadowing terms at nadir
    denom_st0v = np.abs(1.0 + 8.0 * rv0 / (cs * Ftv)) ** 2
    denom_st0h = np.abs(1.0 + 8.0 * rv0 / (cs * Fth)) ** 2
    
    if np.abs(denom_st0v) < 1e-10:
        St0v = 1.0
    else:
        St0v = 1.0 / denom_st0v
    
    if np.abs(denom_st0h) < 1e-10:
        St0h = 1.0
    else:
        St0h = 1.0 / denom_st0h
    
    # Compute series sums
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    temp1 = 1.0
    
    ks_cs = ks * cs
    ks_cs_sq = ks_cs ** 2
    exp_ks_cs_sq = np.exp(-ks_cs_sq)
    
    for n in range(1, min(n_terms, len(spectral_weights)) + 1):
        fn = float(n)
        temp1 *= (1.0 / fn)
        
        a0 = (ks_cs) ** (2.0 * fn)
        weight = spectral_weights[n - 1] if n <= len(spectral_weights) else 0.0
        
        sum1 += temp1 * a0 * weight
        
        # Term for vertical polarization
        term_v = np.abs(Ftv / 2.0 + (2.0 ** (fn + 2.0)) * rv0 / cs * exp_ks_cs_sq) ** 2
        sum2 += temp1 * a0 * term_v * weight
        
        # Term for horizontal polarization
        term_h = np.abs(Fth / 2.0 + (2.0 ** (fn + 2.0)) * rv0 / cs * exp_ks_cs_sq) ** 2
        sum3 += temp1 * a0 * term_h * weight
    
    # Shadowing terms
    if np.abs(sum2) < 1e-10:
        Stv = 0.0
    else:
        Stv = 0.25 * (np.abs(Ftv) ** 2) * sum1 / sum2
    
    if np.abs(sum3) < 1e-10:
        Sth = 0.0
    else:
        Sth = 0.25 * (np.abs(Fth) ** 2) * sum1 / sum3
    
    # Transition functions
    if np.abs(St0v) < 1e-10:
        Tfv = 0.0
    else:
        Tfv = 1.0 - Stv / St0v
    
    if np.abs(St0h) < 1e-10:
        Tfh = 0.0
    else:
        Tfh = 1.0 - Sth / St0h
    
    # Ensure non-negative
    Tfv = max(0.0, float(np.real(Tfv)))
    Tfh = max(0.0, float(np.real(Tfh)))
    
    return Tfv, Tfh
