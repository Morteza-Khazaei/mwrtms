"""I2EM-style transition function for AIEM."""

from __future__ import annotations

import math
import numpy as np

__all__ = ["compute_i2em_transition_function"]


def compute_i2em_transition_function(
    eps_r: complex,
    theta_i: float,
    ks: float,
    cs: float,
    spectral_weights: np.ndarray,
    n_terms: int,
) -> tuple[float, float]:
    """
    Compute I2EM-style transition function for V and H polarizations.
    
    This uses the simpler I2EM transition method which has been shown
    to provide better agreement with NMM3D for co-polarization channels.
    
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
    This implementation follows the I2EM approach which uses a single
    transition factor based on vertical polarization, then applies it
    to both polarizations. This has been empirically shown to reduce
    the systematic bias in co-pol channels.
    
    Key differences from AIEM transition:
    1. Uses single Ft factor (not separate Ftv/Fth)
    2. Simpler computation of St and St0
    3. Same Tf applied to both polarizations
    
    References
    ----------
    Fung, A. K., Li, Z., & Chen, K. S. (1992). Backscattering from a
    randomly rough dielectric surface. IEEE TGRS, 30(2), 356-369.
    """
    si = np.sin(theta_i)
    s2 = si * si
    
    # Nadir reflection coefficients
    # I2EM uses negative rh0 (this is the key difference from AIEM)
    sqrt_er = np.sqrt(eps_r)
    rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
    rh0 = -rv0  # I2EM uses negative for HH
    
    # Transmitted wave vector component
    rt = np.sqrt(eps_r - s2)
    
    # Single transition factor (I2EM style)
    denom = cs * rt
    if np.abs(denom) < 1e-10:
        Ft = 0.0
    else:
        # Note: I2EM uses sin(theta) not sin^2(theta)
        Ft = 8.0 * (rv0 ** 2) * s2 * (cs + rt) / denom
    
    # Compute series sums
    a1 = 0.0
    b1 = 0.0
    
    ks_cs = ks * cs
    ks_cs_sq = ks_cs ** 2
    exp_ks_cs_sq = np.exp(-ks_cs_sq)
    
    for n in range(1, min(n_terms, len(spectral_weights)) + 1):
        fn = float(n)
        weight = spectral_weights[n - 1] if n <= len(spectral_weights) else 0.0
        
        a0 = (ks_cs) ** (2 * fn) / math.factorial(n)
        a1 += a0 * weight
        
        # I2EM uses 2^(n+1) instead of 2^(n+2)
        term = np.abs(Ft / 2.0 + (2.0 ** (fn + 1.0)) * rv0 / cs * exp_ks_cs_sq) ** 2
        b1 += a0 * term * weight
    
    # Shadowing terms (I2EM style)
    if np.abs(b1) < 1e-10:
        St = 0.0
    else:
        St = 0.25 * (np.abs(Ft) ** 2) * a1 / b1
    
    denom_st0 = np.abs(1.0 + 8.0 * rv0 / (cs * Ft)) ** 2
    if np.abs(denom_st0) < 1e-10:
        St0 = 1.0
    else:
        St0 = 1.0 / denom_st0
    
    # Single transition function
    if np.abs(St0) < 1e-10:
        Tf = 0.0
    else:
        Tf = 1.0 - St / St0
    
    # Ensure non-negative
    Tf = max(0.0, float(np.real(Tf)))
    
    # Return same Tf for both polarizations (I2EM approach)
    return Tf, Tf
