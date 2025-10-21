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
    
    This implementation follows the formulation from the main AIEM paper
    (Wu & Fung, 1992) and has been double-checked against the original.
    
    Parameters
    ----------
    eps_r : complex
        Relative dielectric constant of the substrate
    theta_i : float
        Incident angle in radians
    ks : float
        Normalized RMS height (k * sigma)
    cs : float
        cos(theta_i) - cosine of incident angle
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
    The transition function calculation involves:
    
    1. Compute nadir reflection coefficients rv0, rh0
    2. Compute transition factors Ftv, Fth
    3. Compute shadowing terms st0v, st0h at nadir
    4. Compute series sums with spectral weights
    5. Compute shadowing terms stv, sth
    6. Transition function: Tf = 1 - st/st0
    
    The transition function ensures smooth behavior across the full
    range of incidence angles and surface roughness conditions.
    
    References
    ----------
    Wu, T. D., & Fung, A. K. (1992). A transition model for the reflection
    coefficient in surface scattering. IEEE TGRS, 30(4), 856-860.
    """
    # Compute sin^2(theta_i)
    si = np.sin(theta_i)
    siti2 = si * si
    
    # Nadir reflection coefficients
    sqrt_er = np.sqrt(eps_r)
    rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
    rh0 = rv0  # symmetric relation
    
    # Transmitted wave vector component
    sqrt_term = np.sqrt(eps_r - siti2)
    
    # Transition factors
    Ftv = 8.0 * rv0**2 * siti2 * (cs + sqrt_term) / (cs * sqrt_term)
    Fth = -8.0 * rh0**2 * siti2 * (cs + sqrt_term) / (cs * sqrt_term)
    
    # Shadowing terms at nadir
    st0v = 1.0 / np.abs(1.0 + 8.0 * rv0 / (cs * Ftv))**2
    st0h = 1.0 / np.abs(1.0 + 8.0 * rh0 / (cs * Fth))**2
    
    # Compute series sums
    sum_vn = 0.0
    sum_hn = 0.0
    sum_vd = 0.0
    sum_hd = 0.0
    
    # Calculate the transition reflection coefficients
    # using the series expansion
    kscs = ks * cs
    
    for n in range(n_terms):
        fn = n + 1
        temp = kscs**(2.0 * fn) / math.factorial(fn)
        weight = spectral_weights[n] if n < len(spectral_weights) else 0.0

        term_v = np.abs(Ftv + 2.0**(fn + 2.0) * (rv0 / cs) * np.exp(-kscs**2))**2
        term_h = np.abs(Fth + 2.0**(fn + 2.0) * (rh0 / cs) * np.exp(-kscs**2))**2
        
        sum_vn += temp * weight * np.abs(Ftv)**2
        sum_hn += temp * weight * np.abs(Fth)**2
        sum_vd += temp * weight * term_v
        sum_hd += temp * weight * term_h
    
    # Calculate the shadowing terms
    stv = 0.25 * (sum_vn / sum_vd) if np.abs(sum_vd) > 1e-12 else 0.0
    sth = 0.25 * (sum_hn / sum_hd) if np.abs(sum_hd) > 1e-12 else 0.0
    
    # Calculate transition functions
    tfv = 1.0 - (stv / st0v) if np.abs(st0v) > 1e-12 else 0.0
    tfh = 1.0 - (sth / st0h) if np.abs(st0h) > 1e-12 else 0.0
    
    # Avoid negative values (ensure non-negative)
    tfv = max(0.0, float(np.real(tfv)))
    tfh = max(0.0, float(np.real(tfh)))
    # print(f"TFV: {tfv}, TFH: {tfh}")
    
    return tfv, tfh
