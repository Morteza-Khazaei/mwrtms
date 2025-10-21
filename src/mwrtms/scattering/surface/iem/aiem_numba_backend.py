"""Numba-accelerated backend for AIEM multiple scattering computations.

This module provides Numba-jitted kernels for the computationally intensive
parts of the AIEM multiple scattering calculation, including:
- Roughness spectrum computation
- Series summation with factorials
- Coefficient computation (C and B coefficients)
- Integration loops

If Numba is not available, NUMBA_AVAILABLE is False and callers must fall
back to the standard NumPy-based implementation.
"""

from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Fallback decorators
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    def prange(n):
        return range(n)


# ============================================================================
# Numba-Accelerated Roughness Spectrum
# ============================================================================

@njit(cache=True, fastmath=True)
def compute_wn_exponential_numba(u: float, v: float, n: int, sigma2: float, kl: float, two_pi_power: float) -> float:
    """Compute exponential roughness spectrum for order n.
    
    Parameters
    ----------
    u, v : float
        Spectral coordinates
    n : int
        Spectral order
    sigma2 : float
        Variance (σ²)
    kl : float
        Normalized correlation length (k*l)
    two_pi_power : float
        Normalization factor (default 2π for exponential correlation)
        
    Returns
    -------
    float
        Roughness spectrum value
    """
    if n < 1:
        n = 1
    rho_sq = u * u + v * v
    rho = math.sqrt(rho_sq)
    denom = 1.0 + (kl * rho / n) ** 2
    return two_pi_power * sigma2 * (kl / n) ** 2 * denom ** (-1.5)


@njit(cache=True, fastmath=True)
def compute_wn_gaussian_numba(u: float, v: float, n: int, sigma2: float, kl: float) -> float:
    """Compute Gaussian roughness spectrum for order n.
    
    Parameters
    ----------
    u, v : float
        Spectral coordinates
    n : int
        Spectral order
    sigma2 : float
        Variance (σ²)
    kl : float
        Normalized correlation length (k*l)
        
    Returns
    -------
    float
        Roughness spectrum value
    """
    if n < 1:
        n = 1
    factor = kl * kl / n
    exp_arg = -(kl * kl / (4.0 * n)) * (u * u + v * v)
    return sigma2 * (factor / (4.0 * math.pi)) * math.exp(exp_arg)


# ============================================================================
# Numba-Accelerated Series Summation
# ============================================================================

@njit(cache=True, fastmath=True)
def series_sum_exponential_numba(
    coeff_real: float,
    coeff_imag: float,
    arg_x: float,
    arg_y: float,
    sigma2: float,
    kl: float,
    nmax: int,
    two_pi_power: float,
    factorials: np.ndarray
) -> tuple[float, float]:
    """Compute series summation for exponential correlation.
    
    Uses pre-computed factorials and Kahan summation for numerical stability.
    
    Parameters
    ----------
    coeff_real, coeff_imag : float
        Real and imaginary parts of coefficient
    arg_x, arg_y : float
        Spectral arguments
    sigma2 : float
        Surface variance
    kl : float
        Normalized correlation length
    nmax : int
        Maximum order
    two_pi_power : float
        Normalization factor (default 2π for exponential correlation)
    factorials : np.ndarray
        Pre-computed factorials [1!, 2!, ..., nmax!]
        
    Returns
    -------
    tuple[float, float]
        Real and imaginary parts of sum
    """
    result_real = 0.0
    result_imag = 0.0
    c_real = 0.0  # Kahan compensation
    c_imag = 0.0
    
    # Current power of coefficient
    pow_real = coeff_real
    pow_imag = coeff_imag
    
    for n in range(1, nmax + 1):
        # Compute W_n
        wn = compute_wn_exponential_numba(arg_x, arg_y, n, sigma2, kl, two_pi_power)
        
        # Divide by n!
        fact = factorials[n - 1]
        term_real = (pow_real * wn) / fact
        term_imag = (pow_imag * wn) / fact
        
        # Kahan summation for real part
        y_real = term_real - c_real
        t_real = result_real + y_real
        c_real = (t_real - result_real) - y_real
        result_real = t_real
        
        # Kahan summation for imaginary part
        y_imag = term_imag - c_imag
        t_imag = result_imag + y_imag
        c_imag = (t_imag - result_imag) - y_imag
        result_imag = t_imag
        
        # Update power: pow = pow * coeff
        if n < nmax:
            new_real = pow_real * coeff_real - pow_imag * coeff_imag
            new_imag = pow_real * coeff_imag + pow_imag * coeff_real
            pow_real = new_real
            pow_imag = new_imag
    
    return result_real, result_imag


@njit(cache=True, fastmath=True)
def series_sum_gaussian_numba(
    coeff_real: float,
    coeff_imag: float,
    arg_x: float,
    arg_y: float,
    sigma2: float,
    kl: float,
    nmax: int,
    factorials: np.ndarray
) -> tuple[float, float]:
    """Compute series summation for Gaussian correlation.
    
    Parameters
    ----------
    coeff_real, coeff_imag : float
        Real and imaginary parts of coefficient
    arg_x, arg_y : float
        Spectral arguments
    sigma2 : float
        Surface variance
    kl : float
        Normalized correlation length
    nmax : int
        Maximum order
    factorials : np.ndarray
        Pre-computed factorials
        
    Returns
    -------
    tuple[float, float]
        Real and imaginary parts of sum
    """
    result_real = 0.0
    result_imag = 0.0
    
    pow_real = coeff_real
    pow_imag = coeff_imag
    
    for n in range(1, nmax + 1):
        wn = compute_wn_gaussian_numba(arg_x, arg_y, n, sigma2, kl)
        fact = factorials[n - 1]
        
        result_real += (pow_real * wn) / fact
        result_imag += (pow_imag * wn) / fact
        
        if n < nmax:
            new_real = pow_real * coeff_real - pow_imag * coeff_imag
            new_imag = pow_real * coeff_imag + pow_imag * coeff_real
            pow_real = new_real
            pow_imag = new_imag
    
    return result_real, result_imag


# ============================================================================
# Vectorised Series Summation over Precomputed Stacks
# ============================================================================

@njit(parallel=True, fastmath=True)
def series_sum_stack_numba(
    coeff: np.ndarray,
    stack: np.ndarray,
    inv_factorials: np.ndarray,
) -> np.ndarray:
    """Compute series sums for each coefficient using a precomputed W-stack.

    Parameters
    ----------
    coeff : np.ndarray
        Complex coefficient array (...).
    stack : np.ndarray
        Real-valued roughness spectra stack (..., Nmax).
    inv_factorials : np.ndarray
        Array of 1 / n! for n=1..Nmax.

    Returns
    -------
    np.ndarray
        Complex series sum with the same leading shape as ``coeff``.
    """
    leading = coeff.size
    n_terms = stack.shape[-1]
    result = np.empty(leading, dtype=np.complex128)

    coeff_flat = coeff.reshape(leading)
    stack_flat = stack.reshape(leading, n_terms)

    for idx in prange(leading):
        coeff_val = coeff_flat[idx]
        pow_val = coeff_val
        acc = 0.0 + 0.0j
        for n_idx in range(n_terms):
            acc += pow_val * stack_flat[idx, n_idx] * inv_factorials[n_idx]
            pow_val *= coeff_val
        result[idx] = acc

    return result.reshape(coeff.shape)


# ============================================================================
# Numba-Accelerated Integration
# ============================================================================

@njit(cache=True, fastmath=True, parallel=True)
def integrate_2d_weighted_numba(
    integrand_real: np.ndarray,
    integrand_imag: np.ndarray,
    weights_2d: np.ndarray,
    mask: np.ndarray
) -> tuple[float, float]:
    """Perform weighted 2D integration with mask using parallel loops.
    
    Parameters
    ----------
    integrand_real, integrand_imag : np.ndarray
        Real and imaginary parts of integrand (2D arrays)
    weights_2d : np.ndarray
        2D weight array (outer product of 1D weights)
    mask : np.ndarray
        Boolean mask for radiation condition
        
    Returns
    -------
    tuple[float, float]
        Real and imaginary parts of integral
    """
    ni, nj = integrand_real.shape
    total_real = 0.0
    total_imag = 0.0
    
    for i in prange(ni):
        local_real = 0.0
        local_imag = 0.0
        for j in range(nj):
            if mask[i, j]:
                w = weights_2d[i, j]
                local_real += integrand_real[i, j] * w
                local_imag += integrand_imag[i, j] * w
        total_real += local_real
        total_imag += local_imag
    
    return total_real, total_imag


@njit(cache=True, fastmath=True, parallel=True)
def integrate_2d_real_numba(
    integrand: np.ndarray,
    weights_2d: np.ndarray,
    mask: np.ndarray
) -> float:
    """Perform weighted 2D integration for real-valued integrand.
    
    Parameters
    ----------
    integrand : np.ndarray
        Real integrand (2D array)
    weights_2d : np.ndarray
        2D weight array
    mask : np.ndarray
        Boolean mask
        
    Returns
    -------
    float
        Integral value
    """
    ni, nj = integrand.shape
    total = 0.0
    
    for i in prange(ni):
        local_sum = 0.0
        for j in range(nj):
            if mask[i, j]:
                local_sum += integrand[i, j] * weights_2d[i, j]
        total += local_sum
    
    return total


# ============================================================================
# Utility Functions
# ============================================================================

@njit(cache=True, fastmath=True)
def complex_sqrt_numba(real: float, imag: float) -> tuple[float, float]:
    """Compute complex square root.
    
    Parameters
    ----------
    real, imag : float
        Real and imaginary parts of input
        
    Returns
    -------
    tuple[float, float]
        Real and imaginary parts of sqrt
    """
    if imag == 0.0:
        if real >= 0.0:
            return math.sqrt(real), 0.0
        else:
            return 0.0, math.sqrt(-real)
    
    mag = math.sqrt(real * real + imag * imag)
    w = math.sqrt(0.5 * (mag + real))
    
    if real >= 0.0:
        return w, 0.5 * imag / w
    else:
        if imag >= 0.0:
            return 0.5 * imag / w, w
        else:
            return -0.5 * imag / w, -w


@njit(cache=True, fastmath=True)
def safe_divide_complex_numba(
    a_real: float, a_imag: float,
    b_real: float, b_imag: float,
    eps: float = 1e-8
) -> tuple[float, float]:
    """Safe complex division with singularity avoidance.
    
    Parameters
    ----------
    a_real, a_imag : float
        Numerator (real and imaginary)
    b_real, b_imag : float
        Denominator (real and imaginary)
    eps : float
        Safety epsilon
        
    Returns
    -------
    tuple[float, float]
        Result (real and imaginary)
    """
    b_mag_sq = b_real * b_real + b_imag * b_imag
    
    if b_mag_sq < eps * eps:
        b_real += eps
        b_imag += eps
        b_mag_sq = b_real * b_real + b_imag * b_imag
    
    # (a_real + i*a_imag) / (b_real + i*b_imag)
    # = (a_real*b_real + a_imag*b_imag + i*(a_imag*b_real - a_real*b_imag)) / b_mag_sq
    result_real = (a_real * b_real + a_imag * b_imag) / b_mag_sq
    result_imag = (a_imag * b_real - a_real * b_imag) / b_mag_sq
    
    return result_real, result_imag


@njit(cache=True, fastmath=True)
def fresnel_coefficients_numba(
    er_real: float,
    er_imag: float,
    q1_real: float,
    q1_imag: float,
    q2_real: float,
    q2_imag: float
) -> tuple[float, float, float, float]:
    """Compute Fresnel reflection coefficients (Rh and Rv).
    
    Parameters
    ----------
    er_real, er_imag : float
        Relative permittivity
    q1_real, q1_imag : float
        Vertical wavenumber in air
    q2_real, q2_imag : float
        Vertical wavenumber in substrate
        
    Returns
    -------
    tuple[float, float, float, float]
        Rh_real, Rh_imag, Rv_real, Rv_imag
    """
    # Rh = (q1 - q2) / (q1 + q2)
    num_h_real = q1_real - q2_real
    num_h_imag = q1_imag - q2_imag
    den_h_real = q1_real + q2_real
    den_h_imag = q1_imag + q2_imag
    Rh_real, Rh_imag = safe_divide_complex_numba(num_h_real, num_h_imag, den_h_real, den_h_imag)
    
    # Rv = (er*q1 - q2) / (er*q1 + q2)
    er_q1_real = er_real * q1_real - er_imag * q1_imag
    er_q1_imag = er_real * q1_imag + er_imag * q1_real
    num_v_real = er_q1_real - q2_real
    num_v_imag = er_q1_imag - q2_imag
    den_v_real = er_q1_real + q2_real
    den_v_imag = er_q1_imag + q2_imag
    Rv_real, Rv_imag = safe_divide_complex_numba(num_v_real, num_v_imag, den_v_real, den_v_imag)
    
    return Rh_real, Rh_imag, Rv_real, Rv_imag


# ============================================================================
# Pre-computation Utilities
# ============================================================================

def precompute_factorials(nmax: int) -> np.ndarray:
    """Pre-compute factorials up to nmax.
    
    Parameters
    ----------
    nmax : int
        Maximum order
        
    Returns
    -------
    np.ndarray
        Array of factorials [1!, 2!, ..., nmax!]
    """
    factorials = np.zeros(nmax, dtype=np.float64)
    fact = 1.0
    for n in range(1, nmax + 1):
        fact *= n
        factorials[n - 1] = fact
    return factorials


def get_two_pi_power(power: int = 10) -> float:
    """Compute (2π)^power.
    
    Parameters
    ----------
    power : int
        Exponent (default 10)
        
    Returns
    -------
    float
        (2π)^power
    """
    return (2.0 * math.pi) ** power


# ============================================================================
# High-Level Numba Integration Function
# ============================================================================

@njit(cache=True, fastmath=True)
def compute_kirchhoff_integrand_numba(
    U: np.ndarray,
    V: np.ndarray,
    q1: np.ndarray,
    propagator_fp_abs_sq: np.ndarray,
    propagator_fm_abs_sq: np.ndarray,
    propagator_gp_abs_sq: np.ndarray,
    K1: np.ndarray,
    K2: np.ndarray,
    K3: np.ndarray
) -> np.ndarray:
    """Compute Kirchhoff-complementary integrand.
    
    Parameters
    ----------
    U, V : np.ndarray
        Spectral coordinates
    q1 : np.ndarray
        Vertical wavenumber in air
    propagator_*_abs_sq : np.ndarray
        |Propagator|² values
    K1, K2, K3 : np.ndarray
        Kirchhoff-complementary terms
        
    Returns
    -------
    np.ndarray
        Kirchhoff integrand
    """
    ni, nj = U.shape
    result = np.zeros((ni, nj), dtype=np.float64)
    
    for i in range(ni):
        for j in range(nj):
            result[i, j] = (
                propagator_fp_abs_sq[i, j] * K1[i, j] +
                propagator_fm_abs_sq[i, j] * K2[i, j] +
                propagator_gp_abs_sq[i, j] * K3[i, j]
            )
    
    return result


# ============================================================================
# Performance Monitoring
# ============================================================================

def check_numba_performance():
    """Check if Numba is available and working.
    
    Returns
    -------
    dict
        Status information
    """
    status = {
        'numba_available': NUMBA_AVAILABLE,
        'numba_version': None,
        'test_passed': False
    }
    
    if NUMBA_AVAILABLE:
        try:
            import numba
            status['numba_version'] = numba.__version__
            
            # Quick test
            test_result = compute_wn_exponential_numba(1.0, 1.0, 1, 1.0, 1.0, 1.0)
            status['test_passed'] = isinstance(test_result, float)
        except Exception as e:
            status['error'] = str(e)
    
    return status


if __name__ == "__main__":
    # Test the module
    print("AIEM Numba Backend Status:")
    print("-" * 50)
    status = check_numba_performance()
    for key, value in status.items():
        print(f"{key}: {value}")
    
    if NUMBA_AVAILABLE:
        print("\n✅ Numba acceleration is available!")
        print("Expected speedup: 20-100x for multiple scattering")
    else:
        print("\n⚠️  Numba not available - using NumPy fallback")
        print("Install numba: pip install numba")
