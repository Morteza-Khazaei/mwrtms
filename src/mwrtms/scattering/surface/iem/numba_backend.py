"""Numba-accelerated backend for IEM-family computations.

This module provides Numba-jitted kernels for the heavy multiple-scattering
cross-polarization integral in I2EM, plus simple Simpson-rule grids.

If Numba is not available, NUMBA_AVAILABLE is False and the exported
functions should not be used (callers must fall back to the SciPy-based path).
"""
from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover - safe fallback when numba is absent
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        def deco(f):
            return f
        return deco

    def prange(n):  # type: ignore
        return range(n)


from numpy.polynomial.legendre import leggauss


def make_uniform_grid(a: float, b: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Build a Simpson-rule-friendly uniform grid and weights on [a, b].

    Ensures odd number of samples for Simpson's 1/3 rule.
    Returns (points, weights).
    """
    if n % 2 == 0:
        n += 1
    x = np.linspace(a, b, n, dtype=np.float64)
    w = np.ones(n, dtype=np.float64)
    if n >= 3:
        w[1:-1:2] = 4.0
        w[2:-2:2] = 2.0
    w *= (b - a) / (3.0 * (n - 1))
    return x, w


def gauss_legendre_grid(a: float, b: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Build a Gauss–Legendre grid and weights on [a, b]."""
    if n < 2:
        n = 2
    x, w = leggauss(n)
    xm = 0.5 * (b - a)
    xc = 0.5 * (b + a)
    pts = xm * x + xc
    wts = xm * w
    return pts.astype(np.float64), wts.astype(np.float64)


@njit(cache=True)
def cmath_sqrt(z: complex) -> complex:
    # Numba-friendly complex sqrt
    return z ** 0.5


@njit(cache=True, fastmath=False)
def _xpol_vhmnsum(ks2: float, cs2: float, kl2: float, rx: float, ry: float, s: float, n_spec: int, fact: np.ndarray) -> float:
    """Compute the double spectral sum with Kahan-compensated accumulation."""
    v = 0.0
    c = 0.0
    for m in range(1, n_spec + 1):
        for n in range(1, n_spec + 1):
            den_wn = (n * n) + kl2 * ((rx - s) * (rx - s) + ry * ry)
            wn = n * kl2 / (den_wn * math.sqrt(den_wn))
            den_wm = (m * m) + kl2 * ((rx + s) * (rx + s) + ry * ry)
            wm = m * kl2 / (den_wm * math.sqrt(den_wm))
            p = (ks2 * cs2) ** (m + n)
            term = p * (wn * wm) / (fact[n - 1] * fact[m - 1])
            y = term - c
            t = v + y
            c = (t - v) - y
            v = t
    return v


@njit(cache=True, fastmath=True, parallel=True)
def xpol_integrate_numba(
    ks: float,
    kl: float,
    theta: float,
    er: complex,
    rss: float,
    n_spec: int,
    r_grid: np.ndarray,
    r_w: np.ndarray,
    phi_grid: np.ndarray,
    phi_w: np.ndarray,
) -> float:
    """Compute the cross-pol integral using Numba over precomputed grids.

    This implements the same integrand as _xpol_integral in i2em.py, including
    the 1e5 scaling inside the integrand, but integrates via Simpson weights.
    The caller should apply the same backscatter shadowing factor times 1e-5,
    which cancels the scale and mirrors the original code path exactly.
    """
    cs = math.cos(theta)
    if abs(cs) < 1e-6:
        cs = 1e-6
    s = math.sin(theta + 0.001)  # restore small offset for stability
    ks2 = ks * ks
    kl2 = kl * kl
    cs2 = cs * cs

    # factorials 1..n_spec
    fact = np.empty(n_spec, dtype=np.float64)
    acc = 1.0
    for i in range(n_spec):
        acc *= (i + 1)
        fact[i] = acc

    total = 0.0

    for ip in prange(phi_grid.shape[0]):
        phi = float(phi_grid[ip])
        wphi = float(phi_w[ip])
        sp = math.sin(phi)
        cp = math.cos(phi)
        for ir in range(r_grid.shape[0]):
            r = float(r_grid[ir])
            wr = float(r_w[ir])

            # Reflection coefficient difference and field coefficients
            rt = cmath_sqrt(er - s * s)
            rv = (er * cs - rt) / (er * cs + rt)
            rh = (cs - rt) / (cs + rt)
            rvh = (rv - rh) * 0.5
            rp = 1.0 + rvh
            rm = 1.0 - rvh
            t = 1.0001 - r * r
            q = math.sqrt(t) if t > 0.0 else 0.0
            qt = cmath_sqrt(er - r * r)
            a = rp / q
            b = rm / q
            c = rp / qt
            d = rm / qt
            B3 = (r * sp) * (r * cp) / cs
            fvh1 = (b - c) * (1.0 - 3.0 * rvh) - (b - c / er) * rp
            fvh2 = (a - d) * (1.0 + 3.0 * rvh) - (a - d * er) * rm
            Fvh = abs((fvh1 + fvh2) * B3) ** 2

            # Shadowing for multiple scattering
            au = q / (r * 1.4142135623730951 * rss)
            fsh = (0.2821 / au) * math.exp(-au * au) - 0.5 * math.erfc(au)
            sha = 1.0 / (1.0 + fsh)

            # Spectral terms and double sum
            rx = r * cp
            ry = r * sp
            vhmn = _xpol_vhmnsum(ks2, cs2, kl2, rx, ry, s, n_spec, fact)

            acc_factor = math.exp(-2.0 * ks2 * cs2) / (16.0 * math.pi)
            VH = 4.0 * acc_factor * Fvh * vhmn * r

            # Match the original integrand's 1e5 scaling
            total += (VH * sha * 1e5) * wr * wphi

    return total
