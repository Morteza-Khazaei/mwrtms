"""AIEM-specific roughness spectrum computation."""

from __future__ import annotations

import numpy as np
from scipy.special import gamma, kv

__all__ = ["compute_aiem_spectrum"]


def compute_aiem_spectrum(
    kl: float,
    K: float,
    n: int,
    correlation_type: str,
    power_exponent: float = 1.5,
    kl_x: float | None = None,
    kl_y: float | None = None,
    kx: float | None = None,
    ky: float | None = None,
) -> float:
    """
    Compute AIEM roughness spectrum for order n.
    
    Parameters
    ----------
    kl : float
        Normalized correlation length (k * L) for isotropic case
    K : float
        Spatial frequency magnitude
    n : int
        Spectral order (1, 2, 3, ...)
    correlation_type : str
        Surface correlation function type:
        - 'gaussian' or '1': Gaussian correlation
        - 'exponential' or '2': Exponential correlation
        - 'powerlaw' or '3': 1.5-power correlation
    power_exponent : float, optional
        Exponent for power-law correlation (default: 1.5)
    kl_x : float, optional
        Normalized correlation length along x-axis (k * Lx) for anisotropic case
    kl_y : float, optional
        Normalized correlation length along y-axis (k * Ly) for anisotropic case
    kx : float, optional
        Spatial frequency x-component (required for anisotropic case)
    ky : float, optional
        Spatial frequency y-component (required for anisotropic case)
    
    Returns
    -------
    spectrum : float
        Roughness spectrum value W_n(K)
    
    Notes
    -----
    From AIEM.m and Yang & Chen (2019):
    
    Gaussian (itype=1):
        W_n = (kl^2 / (2*n)) * exp(-K^2 / (4*n))
    
    Exponential (itype=2):
        Isotropic: W_n = (kl/n)^2 * (1 + (K/n)^2)^(-1.5)
        Anisotropic: W_n(K,φ) = (kl(φ)/n)^2 * (1 + (K*kl(φ)/n)^2)^(-1.5)
                     where kl(φ) = kl_x * cos²(φ) + kl_y * sin²(φ)
    
    1.5-power (itype=3):
        W_n = (kl^2 * (K/2)^e * K_m(K)) / gamma(y)
        where e = 1.5*n - 1, y = 1.5*n, m = 1.5*n - 1
        K_m is the modified Bessel function of the second kind
    
    References
    ----------
    Chen, K. S., et al. (2003). Emission of rough surfaces calculated by the
    integral equation method with comparison to three-dimensional moment method
    simulations. IEEE TGRS, 41(1), 90-101.
    
    Yang, Y., & Chen, K. S. (2019). Polarized backscattering from spatially
    anisotropic rough surface. IEEE TGRS, 57(9), 6608-6618.
    """
    fn = float(n)
    
    # Normalize correlation type
    corr_type = str(correlation_type).lower()
    if corr_type in ('1', 'gaussian', 'gauss'):
        # Gaussian correlation function (isotropic only for now)
        kl2 = kl * kl
        spectrum = (kl2 / (2.0 * fn)) * np.exp(-(K * K) / (4.0 * fn))
    
    elif corr_type in ('2', 'exponential', 'exp'):
        # Exponential correlation function
        # Check if anisotropic parameters are provided
        is_anisotropic = (
            kl_x is not None
            and kl_y is not None
            and kx is not None
            and ky is not None
            and not np.isclose(kl_x, kl_y)
        )
        
        if is_anisotropic:
            # Anisotropic case: kl(φ) = kl_x * cos²(φ) + kl_y * sin²(φ)
            if K > 1e-12:
                cos_phi = kx / K
                sin_phi = ky / K
                kl_phi = kl_x * cos_phi**2 + kl_y * sin_phi**2
            else:
                # At K=0, use average
                kl_phi = 0.5 * (kl_x + kl_y)
            spectrum = (kl_phi / fn) ** 2 * (1.0 + (K / fn) ** 2) ** (-1.5)
        else:
            # Isotropic case
            spectrum = (kl / fn) ** 2 * (1.0 + (K / fn) ** 2) ** (-1.5)
    
    elif corr_type in ('3', 'powerlaw', 'power', '1.5-power', 'xpower'):
        # 1.5-power (transformed exponential) correlation function
        if np.isclose(K, 0.0):
            # Special case: K = 0
            spectrum = kl2 / (3.0 * fn - 2.0)
        else:
            e = power_exponent * fn - 1.0
            y = power_exponent * fn
            m = power_exponent * fn - 1.0
            
            # Use log-domain to avoid overflow
            try:
                log_gamma = np.log(gamma(y))
                log_bessel = np.log(kv(m, K))
                log_out = np.log(kl2) + e * np.log(K / 2.0)
                spectrum = np.exp(log_out + log_bessel - log_gamma)
            except (ValueError, RuntimeWarning):
                # Fallback for numerical issues
                spectrum = kl2 * (K / 2.0) ** e * kv(m, K) / gamma(y)
    
    else:
        raise ValueError(
            f"Unknown correlation type '{correlation_type}'. "
            f"Supported: 'gaussian' (1), 'exponential' (2), 'powerlaw' (3)"
        )
    
    return float(spectrum)
