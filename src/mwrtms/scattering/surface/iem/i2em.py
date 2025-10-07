"""I2EM backscatter implementation built on mwRTMs utilities."""

from __future__ import annotations

import math

from typing import Optional

import numpy as np
from scipy.special import erfc
from scipy.integrate import dblquad
from .base import IEMBase, SurfaceRoughnessParameters
from ....core import PolarizationState
from ....medium import Medium

__all__ = ["I2EMModel"]


class I2EMModel(IEMBase):
    """Python translation of the single-scale I2EM backscatter model."""

    MODEL_NAME = "I2EM"

    def __init__(
        self,
        wave,
        geometry,
        surface,
        *,
        correlation_type: str = "exponential",
        correlation_length_m: Optional[float] = None,
        spectral_terms: Optional[int] = None,
        power_exponent: float = 1.5,
        auto_terms: bool = True,
    ) -> None:
        super().__init__(
            wave,
            geometry,
            surface,
            correlation_type=correlation_type,
            correlation_length_m=correlation_length_m,
            spectral_terms=spectral_terms,
            power_exponent=power_exponent,
            auto_terms=auto_terms,
        )
        # Numba backend toggles and grid density (Simpson-friendly odd sizes)
        # Phase 0: force SciPy dblquad for HV by disabling Numba path
        self._use_numba = True
        self._numba_r_points = 81
        self._numba_phi_points = 97

    def _compute_channel(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
        params: SurfaceRoughnessParameters,
    ) -> float:
        if polarization.is_copol:
            return self._compute_copol_channel(medium_above, medium_below, polarization, params)
        if polarization.is_crosspol:
            return self._compute_crosspol_channel(medium_above, medium_below, polarization, params)

        return 0.0

    def _compute_crosspol_channel(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
        params: SurfaceRoughnessParameters,
    ) -> float:
        """Computes the cross-polarization channel (HV or VH)."""
        eps2 = medium_below.permittivity(self._wave.frequency_hz)

        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_s_rad

        ks = params.ks
        kl = params.kl

        # For backscatter, theta_s = theta_i, phi_s = pi
        # The integral is over the surface plane, so we use dummy vars r, phi
        # The MATLAB reference integrates over a normalized slope domain.
        integral = None
        if getattr(self, "_use_numba", True):
            try:
                from .numba_backend import NUMBA_AVAILABLE, make_uniform_grid, gauss_legendre_grid, xpol_integrate_numba_strict
                if NUMBA_AVAILABLE:
                    r_grid, r_w = make_uniform_grid(1e-4, 1.0, getattr(self, "_numba_r_points", 81))
                    phi_grid, phi_w = make_uniform_grid(0.0, np.pi, getattr(self, "_numba_phi_points", 97))
                    integral = xpol_integrate_numba_strict(
                        ks, kl, theta_i, complex(eps2), params.rms_slope,
                        len(params.spectral_weights),
                        r_grid, r_w, phi_grid, phi_w,
                    )
            except Exception:
                integral = None

        if integral is None:
            # Fallback: SciPy dblquad with the original Python integrand
            integral, _ = dblquad(
                _xpol_integral,
                0, np.pi,          # phi integration limits
                lambda phi: 1e-4, 1, # r integration limits
                args=(ks, kl, theta_i, eps2, self._correlation_type, params.rms_slope, len(params.spectral_weights))
            )

        # The integrand already includes multiple-scattering shadowing; only un-scale.
        sigma_vh = 1e-5 * integral

        # Under reciprocity for backscatter, sigma_hv = sigma_vh
        return float(sigma_vh)

    def _compute_copol_channel(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
        params: SurfaceRoughnessParameters,
    ) -> float:
        """Computes the co-polarization channel (VV or HH)."""
        eps2 = medium_below.permittivity(self._wave.frequency_hz)
        k = self._wave.wavenumber

        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_s_rad
        phi_s = self._geometry.phi_s_rad
        phi_i = 0.0

        # Numerically stable trigonometry without angular offsets
        cs = np.cos(theta_i)
        if np.abs(cs) < 1e-6:
            cs = 1e-6
        s = np.sin(theta_i)
        s2 = s ** 2

        ss = np.sin(theta_s)
        css = np.cos(theta_s)
        cf = np.cos(phi_i)
        cfs = np.cos(phi_s)
        sfs = np.sin(phi_s)

        rt = np.sqrt(eps2 - s2)
        Rvi = (eps2 * cs - rt) / (eps2 * cs + rt) if (eps2 * cs + rt) != 0 else np.inf
        Rhi = (cs - rt) / (cs + rt) if (cs + rt) != 0 else np.inf
        
        rv0 = (np.sqrt(eps2) - 1.0) / (np.sqrt(eps2) + 1.0) if (np.sqrt(eps2) + 1.0) != 0 else np.inf
        rh0 = -rv0

        ks = params.ks
        kz = k * cs
        ksz = k * css

        weights = params.spectral_weights
        ss_vec = params.rms_slope
        n_terms = len(weights)

        Ft = 8.0 * rv0 ** 2 * ss * (cs + np.sqrt(eps2 - s2)) / (cs * np.sqrt(eps2 - s2)) if (cs * np.sqrt(eps2 - s2)) != 0 else np.inf
        a1 = 0.0
        b1 = 0.0
        for n, weight in enumerate(weights, start=1):
            a0 = (ks * cs) ** (2 * n) / math.factorial(n)
            a1 += a0 * weight
            term = np.abs(Ft / 2.0 + 2.0 ** (n + 1) * rv0 / cs * np.exp(-(ks * cs) ** 2)) ** 2
            b1 += a0 * term * weight

        St = 0.25 * (np.abs(Ft) ** 2) * a1 / b1
        St0 = 1.0 / (np.abs(1.0 + 8.0 * rv0 / (cs * Ft)) ** 2) if (cs * Ft) != 0 else np.inf
        Tf = 1.0 - St / St0 if St0 != 0 else 0.0

        if np.isclose(theta_i, theta_s) and np.isclose((phi_s % (2 * np.pi)), np.pi):
            # Backscatter case
            Rvt = Rvi + (rv0 - Rvi) * Tf
            Rht = Rhi + (rh0 - Rhi) * Tf
        else:
            # Bistatic case: use averaged reflectivities
            sigx = 1.1 * params.sigma_m / params.correlation_length_m
            sigy = sigx
            xxx = 3 * sigx
            
            rav_res, _ = dblquad(lambda Zy, Zx: _rav_integration(Zx, Zy, cs, s, eps2, s2, sigx, sigy), -xxx, xxx, lambda x: -xxx, lambda x: xxx)
            rah_res, _ = dblquad(lambda Zy, Zx: _rah_integration(Zx, Zy, cs, s, eps2, s2, sigx, sigy), -xxx, xxx, lambda x: -xxx, lambda x: xxx)

            Rvt = rav_res / (2 * np.pi * sigx * sigy)
            Rht = rah_res / (2 * np.pi * sigx * sigy)
        
        fvv = 2.0 * Rvt * (s * ss - (1.0 + cs * css) * cfs) / (cs + css) if (cs + css) != 0 else np.inf
        fhh = -2.0 * Rht * (s * ss - (1.0 + cs * css) * cfs) / (cs + css) if (cs + css) != 0 else np.inf

        Fvvupi, Fhhupi = _fppupdn_is(+1, 1, Rvi, Rhi, eps2, k, kz, ksz, s, cs, ss, css, cf, cfs, sfs)
        Fvvups, Fhhups = _fppupdn_is(+1, 2, Rvi, Rhi, eps2, k, kz, ksz, s, cs, ss, css, cf, cfs, sfs)
        Fvvdni, Fhhdni = _fppupdn_is(-1, 1, Rvi, Rhi, eps2, k, kz, ksz, s, cs, ss, css, cf, cfs, sfs)
        Fvvdns, Fhhdns = _fppupdn_is(-1, 2, Rvi, Rhi, eps2, k, kz, ksz, s, cs, ss, css, cf, cfs, sfs)
        
        qi = k * cs
        qs = k * css

        Ivv = np.zeros(n_terms, dtype=np.complex128)
        Ihh = np.zeros(n_terms, dtype=np.complex128)
        sigma_m_sq = params.sigma_m ** 2

        for n, weight in enumerate(weights, start=1):
            idx = n - 1
            exp_idx = max(0, n - 1)

            term_vv = (
                0.25 * (
                    Fvvupi * (ksz - qi) ** exp_idx * np.exp(-sigma_m_sq * (qi**2 - qi * (ksz - kz))) +
                    Fvvdni * (ksz + qi) ** exp_idx * np.exp(-sigma_m_sq * (qi**2 + qi * (ksz - kz))) +
                    Fvvups * (kz + qs) ** exp_idx * np.exp(-sigma_m_sq * (qs**2 - qs * (ksz - kz))) +
                    Fvvdns * (kz - qs) ** exp_idx * np.exp(-sigma_m_sq * (qs**2 + qs * (ksz - kz)))
                )
            )
            Ivv[idx] = (kz + ksz) ** n * fvv * np.exp(-sigma_m_sq * kz * ksz) + term_vv

            term_hh = (
                0.25 * (
                    Fhhupi * (ksz - qi) ** exp_idx * np.exp(-sigma_m_sq * (qi**2 - qi * (ksz - kz))) +
                    Fhhdni * (ksz + qi) ** exp_idx * np.exp(-sigma_m_sq * (qi**2 + qi * (ksz - kz))) +
                    Fhhups * (kz + qs) ** exp_idx * np.exp(-sigma_m_sq * (qs**2 - qs * (ksz - kz))) +
                    Fhhdns * (kz - qs) ** exp_idx * np.exp(-sigma_m_sq * (qs**2 + qs * (ksz - kz)))
                )
            )
            Ihh[idx] = (kz + ksz) ** n * fhh * np.exp(-sigma_m_sq * kz * ksz) + term_hh

        is_back = (np.isclose(theta_i, theta_s) and np.isclose((phi_s % (2 * np.pi)), np.pi))
        shdw = _shadowing_factor(theta_i, theta_s, ss_vec) if is_back else 1.0

        sigma_vv = 0.0
        sigma_hh = 0.0
        for n, weight in enumerate(weights, start=1):
            idx = n - 1
            a0 = weight / math.factorial(n) * (params.sigma_m ** (2 * n))
            sigma_vv += np.abs(Ivv[idx]) ** 2 * a0
            sigma_hh += np.abs(Ihh[idx]) ** 2 * a0

        factor = shdw * k**2 / 2.0 * np.exp(-sigma_m_sq * (kz**2 + ksz**2))
        sigma_vv *= factor
        sigma_hh *= factor

        return float(sigma_vv if polarization == PolarizationState.VV else sigma_hh)



def _shadowing_factor(theta_i: float, theta_s: float, rms_slope: float) -> float:
    # Backscatter shadowing factor computed for arbitrary incidence angles.
    # Caller should apply only in true backscatter geometry.
    tan_ti = np.tan(theta_i)
    tan_ts = np.tan(theta_s)
    if np.abs(tan_ti) < 1e-6 or np.abs(tan_ts) < 1e-6 or rms_slope <= 0:
        return 1.0

    ct = 1.0 / tan_ti
    cts = 1.0 / tan_ts

    ctorslp = ct / (np.sqrt(2.0) * rms_slope)
    ctsorslp = cts / (np.sqrt(2.0) * rms_slope)

    shadf = 0.5 * (np.exp(-ctorslp**2) / (np.sqrt(np.pi) * ctorslp) - erfc(ctorslp))
    shadfs = 0.5 * (np.exp(-ctsorslp**2) / (np.sqrt(np.pi) * ctsorslp) - erfc(ctsorslp))

    denom = 1.0 + shadf + shadfs
    return 1.0 / denom if denom != 0.0 else 1.0


def _fppupdn_is(
    ud: int,
    is_case: int,
    Rvi,
    Rhi,
    er,
    k,
    kz,
    ksz,
    s,
    cs,
    ss,
    css,
    cf,
    cfs,
    sfs,
):
    k2 = k**2
    term1 = s * cf
    term2 = ss * cfs
    term3 = ss * sfs

    if is_case == 1:
        Gq = ud * kz
        Gqt = ud * k * np.sqrt(er - s**2)
        q = ud * kz

        c11 = k * cfs * (ksz - q)
        c21 = cs * (cfs * (k2 * term1 * (term2 - term1) + Gq * (k * css - q)) + k2 * cf * s * ss * sfs**2)
        c31 = k * s * (term1 * cfs * (k * css - q) - Gq * (cfs * (term2 - term1) + ss * sfs**2))
        c41 = k * cs * (cfs * css * (k * css - q) + k * ss * (term2 - term1))
        c51 = Gq * (cfs * css * (q - k * css) - k * ss * (term2 - term1))

        c12 = k * cfs * (ksz - q)
        c22 = cs * (cfs * (k2 * term1 * (term2 - term1) + Gqt * (k * css - q)) + k2 * cf * s * ss * sfs**2)
        c32 = k * s * (term1 * cfs * (k * css - q) - Gqt * (cfs * (term2 - term1) - ss * sfs**2))
        c42 = k * cs * (cfs * css * (k * css - q) + k * ss * (term2 - term1))
        c52 = Gqt * (cfs * css * (q - k * css) - k * ss * (term2 - term1))
    else:
        Gq = ud * ksz
        Gqt = ud * k * np.sqrt(er - ss**2)
        q = ud * ksz

        c11 = k * cfs * (kz + q)
        c21 = Gq * (cfs * (cs * (kz + q) - k * s * (term2 - term1)) - k * s * ss * sfs**2)
        c31 = k * ss * (k * cs * (term2 - term1) + s * (kz + q))
        c41 = k * css * (cfs * (cs * (kz + q) - k * s * (term2 - term1)) - k * s * ss * sfs**2)
        c51 = -css * (k2 * ss * (term2 - term1) + Gq * cfs * (kz + q))

        c12 = k * cfs * (kz + q)
        c22 = Gqt * (cfs * (cs * (kz + q) - k * s * (term2 - term1)) - k * s * ss * sfs**2)
        c32 = k * ss * (k * cs * (term2 - term1) + s * (kz + q))
        c42 = k * css * (cfs * (cs * (kz + q) - k * s * (term2 - term1)) - k * s * ss * sfs**2)
        c52 = -css * (k2 * ss * (term2 - term1) + Gqt * cfs * (kz + q))

    qz = kz
    qt = k * np.sqrt(er - s ** 2)

    vv = (
        (1 + Rvi) * (-(1 - Rvi) * c11 / qz + (1 + Rvi) * c12 / qt) +
        (1 - Rvi) * ((1 - Rvi) * c21 / qz - (1 + Rvi) * c22 / qt) +
        (1 + Rvi) * ((1 - Rvi) * c31 / qz - (1 + Rvi) * c32 / (er * qt)) +
        (1 - Rvi) * ((1 + Rvi) * c41 / qz - er * (1 - Rvi) * c42 / qt) +
        (1 + Rvi) * ((1 + Rvi) * c51 / qz - (1 - Rvi) * c52 / qt)
    )

    hh = (
        (1 + Rhi) * ((1 - Rhi) * c11 / qz - er * (1 + Rhi) * c12 / qt) -
        (1 - Rhi) * ((1 - Rhi) * c21 / qz - (1 + Rhi) * c22 / qt) -
        (1 + Rhi) * ((1 - Rhi) * c31 / qz - (1 + Rhi) * c32 / qt) -
        (1 - Rhi) * ((1 + Rhi) * c41 / qz - (1 - Rhi) * c42 / qt) -
        (1 + Rhi) * ((1 + Rhi) * c51 / qz - (1 - Rhi) * c52 / qt)
    )

    return vv, hh


def _rav_integration(Zx, Zy, cs, s, er, s2, sigx, sigy):
    A = cs + Zx * s
    B = er * (1 + Zx**2 + Zy**2)
    CC = s2 - 2 * Zx * s * cs + Zx**2 * cs**2 + Zy**2
    Rv = (er * A - np.sqrt(B - CC)) / (er * A + np.sqrt(B - CC))
    pd = np.exp(-Zx**2 / (2 * sigx**2) - Zy**2 / (2 * sigy**2))
    return Rv * pd


def _rah_integration(Zx, Zy, cs, s, er, s2, sigx, sigy):
    A = cs + Zx * s
    B = er * (1 + Zx**2 + Zy**2)
    CC = s2 - 2 * Zx * s * cs + Zx**2 * cs**2 + Zy**2
    Rh = (A - np.sqrt(B - CC)) / (A + np.sqrt(B - CC))
    pd = np.exp(-Zx**2 / (2 * sigx**2) - Zy**2 / (2 * sigy**2))
    return Rh * pd


def _spectrm1(corr_type: str, r, kl, n_terms=10):
    """First-order surface spectrum, adapted for IEMX_model.m logic."""
    # This function is not directly used by the new _xpol_integral_func_matlab
    # but is kept for potential future use or other IEM variants.
    if corr_type == "exponential":
        return kl * (kl**2 + (r / kl) ** 2) ** -1.5
    elif corr_type == "gaussian":
        return (kl**2 / 2) * np.exp(-((r / 2) ** 2))
    else:
        return kl * (kl**2 + (r / kl) ** 2) ** -1.5


def _spectrm2(corr_type: str, r, kl, n_terms=10):
    """Second-order surface spectrum, adapted for IEMX_model.m logic."""
    # This function is not directly used by the new _xpol_integral_func_matlab
    if corr_type == "exponential":
        return (kl**2) * (kl**2 + (r / kl) ** 2) ** -1.5
    elif corr_type == "gaussian":
        return (kl**2 / 2) * np.exp(-((r / 2) ** 2)) * (2 - (r / kl) ** 2)
    else:
        return (kl**2) * (kl**2 + (r / kl) ** 2) ** -1.5

def _xpol_integral(r, phi, ks, kl, theta, er, corr_type, rss, n_spec):
    """
    Integrand for cross-polarization, translated from IEMX_model.m.
    Note: 'r' and 'phi' here are normalized slope variables, not polar coordinates.
    """
    cs = np.cos(theta)
    if np.abs(cs) < 1e-6:
        cs = 1e-6
    s = np.sin(theta + 0.001)
    ks2 = ks**2
    kl2 = kl**2
    cs2 = cs**2

    # Reflection coefficient difference
    rt = np.sqrt(er - s**2)
    rv = (er * cs - rt) / (er * cs + rt)
    rh = (cs - rt) / (cs + rt)
    rvh = (rv - rh) / 2.0

    # Field coefficients
    rp = 1 + rvh
    rm = 1 - rvh
    t = 1.0001 - r**2
    q = np.sqrt(t) if t > 0.0 else 0.0
    qt = np.sqrt(er - r**2)
    a = rp / q
    b = rm / q
    c = rp / qt
    d = rm / qt
    B3 = r * np.sin(phi) * r * np.cos(phi) / cs
    fvh1 = (b - c) * (1 - 3 * rvh) - (b - c / er) * rp
    fvh2 = (a - d) * (1 + 3 * rvh) - (a - d * er) * rm
    Fvh = np.abs((fvh1 + fvh2) * B3)**2

    # Shadowing for multiple scattering
    au = q / r / 1.414 / rss
    fsh = (0.2821 / au) * np.exp(-au**2) - 0.5 * erfc(au)
    sha = 1.0 / (1.0 + fsh)

    # Surface spectra Wn and Wm
    rx = r * np.cos(phi)
    ry = r * np.sin(phi)
    n_range = np.arange(1, n_spec + 1)
    
    # Wn term
    den_wn = n_range[:, np.newaxis]**2 + kl2 * ((rx - s)**2 + ry**2)
    wn = n_range[:, np.newaxis] * kl2 / (den_wn**1.5)
    
    # Wm term
    den_wm = n_range[:, np.newaxis]**2 + kl2 * ((rx + s)**2 + ry**2)
    wm = n_range[:, np.newaxis] * kl2 / (den_wm**1.5)

    # Double summation
    factorials = np.array([math.factorial(n) for n in n_range])
    vhmnsum = np.sum(
        (ks2 * cs2)**(n_range[:, np.newaxis] + n_range[np.newaxis, :]) *
        (wn[:, :, np.newaxis] * wm[:, np.newaxis, :]).T /
        (factorials[:, np.newaxis] * factorials[np.newaxis, :])
    )

    # Final integrand value
    acc = np.exp(-2 * ks2 * cs2) / (16 * np.pi)
    VH = 4 * acc * Fvh * vhmnsum * r
    
    # Rescale to help dblquad convergence, as in MATLAB
    return VH * sha * 1e5
