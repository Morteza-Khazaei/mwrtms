# cython: boundscheck=False
# cython: wraparound=False
# cython: initializedcheck=False
# cython: nonecheck=False
# cython: cdivision=True

"""Cython implementation of the I2EM bistatic scattering model."""

import math
import numpy as np

cimport cython
cimport numpy as np

from libc.math cimport tgamma
from scipy.special import erfc, gamma, kv as besselk, jv as besselj
from scipy.integrate import quad

ctypedef np.float64_t float64_t
ctypedef np.complex128_t complex128_t


def _roughness_spectrum_case4(double wvnb, double L, int Ts, double xx):
    wn = np.zeros(Ts, dtype=np.float64)
    for n in range(1, Ts + 1):
        def integrand(z, n=n):
            return math.exp(-abs(z) ** xx) * besselj(0, z * wvnb * L / n ** (1.0 / xx)) * z
        integral, _ = quad(integrand, 0.0, 9.0)
        wn[n - 1] = L * L / n ** (2.0 / xx) * integral
    return wn


def _compute_average_reflectivities(double cs, double s, complex er, double s2,
                                    double sigx, double sigy, double xxx):
    def _rav(Zx, Zy):
        return Rav_integration(Zx, Zy, cs, s, er, s2, sigx, sigy)

    def _rah(Zx, Zy):
        return Rah_integration(Zx, Zy, cs, s, er, s2, sigx, sigy)

    def _integrate(func):
        def real_inner(Zx, f=func):
            val, _ = quad(lambda Zy: float(f(Zx, Zy).real), -xxx, xxx)
            return val

        def imag_inner(Zx, f=func):
            val, _ = quad(lambda Zy: float(f(Zx, Zy).imag), -xxx, xxx)
            return val

        real_val, _ = quad(real_inner, -xxx, xxx)
        imag_val, _ = quad(imag_inner, -xxx, xxx)
        return real_val + 1j * imag_val

    factor = 1.0 / (2.0 * math.pi * sigx * sigy)
    Rav = _integrate(_rav) * factor
    Rah = _integrate(_rah) * factor
    return Rav, Rah

cdef inline double factorial_double(int n):
    """Return factorial(n) as double using Gamma function."""
    return tgamma(n + 1.0)

cdef tuple roughness_spectrum(int sp, double xx, double wvnb, double sig,
                              double L, int Ts):
    cdef np.ndarray[np.float64_t, ndim=1] wn = np.zeros(Ts, dtype=np.float64)
    cdef int n
    cdef double rss
    cdef double arg, nu

    if sp == 1:
        for n in range(1, Ts + 1):
            wn[n - 1] = (L * L) / (n * n) * (1.0 + (wvnb * L / n) ** 2) ** (-1.5)
        rss = sig / L
    elif sp == 2:
        for n in range(1, Ts + 1):
            wn[n - 1] = L * L / (2.0 * n) * math.exp(-(wvnb * L) ** 2 / (4.0 * n))
        rss = math.sqrt(2.0) * sig / L
    elif sp == 3:
        for n in range(1, Ts + 1):
            if wvnb == 0.0:
                wn[n - 1] = (L * L) / (3.0 * n - 2.0)
            else:
                arg = wvnb * L
                nu = 1.0 - xx * n
                wn[n - 1] = (
                    (L * L) * arg ** (-nu) * besselk(nu, arg) /
                    (2.0 ** (xx * n - 1.0) * gamma(xx * n))
                )
        rss = math.sqrt(2.0 * xx) * sig / L
    elif sp == 4:
        wn = np.asarray(_roughness_spectrum_case4(wvnb, L, Ts, xx), dtype=np.float64)
        rss = sig / L
    else:
        raise ValueError("Unsupported sp value, must be 1-4")

    return wn, rss

cdef complex Rav_integration(double Zx, double Zy, double cs, double s,
                             complex er, double s2, double sigx, double sigy):
    A = cs + Zx * s
    B = er * (1.0 + Zx * Zx + Zy * Zy)
    CC = s2 - 2.0 * Zx * s * cs + Zx * Zx * cs * cs + Zy * Zy
    root_term = np.sqrt(B - CC)
    Rv = (er * A - root_term) / (er * A + root_term)
    pd = math.exp(-Zx * Zx / (2.0 * sigx * sigx) - Zy * Zy / (2.0 * sigy * sigy))
    return Rv * pd

cdef complex Rah_integration(double Zx, double Zy, double cs, double s,
                             complex er, double s2, double sigx, double sigy):
    A = cs + Zx * s
    B = er * (1.0 + Zx * Zx + Zy * Zy)
    CC = s2 - 2.0 * Zx * s * cs + Zx * Zx * cs * cs + Zy * Zy
    root_term = np.sqrt(B - CC)
    Rh = (A - root_term) / (A + root_term)
    pd = math.exp(-Zx * Zx / (2.0 * sigx * sigx) - Zy * Zy / (2.0 * sigy * sigy))
    return Rh * pd

cdef tuple Fppupdn_is_calculations(int ud, int is_mode, complex Rvi, complex Rhi,
                                   complex er, double k, double kz, double ksz,
                                   double s, double cs, double ss, double css,
                                   double cf, double cfs, double sfs):
    if is_mode == 1:
        Gq = ud * kz
        Gqt = ud * k * np.sqrt(er - s * s)
        q = ud * kz
    elif is_mode == 2:
        Gq = ud * ksz
        Gqt = ud * k * np.sqrt(er - ss * ss)
        q = ud * ksz
    else:
        raise ValueError("is_mode must be 1 or 2")

    k2 = k * k
    term1 = s * cf
    term2 = ss * cfs

    c11 = k * cfs * (ksz - q)
    c12 = k * cfs * (ksz - q)

    if is_mode == 1:
        c21 = cs * (cfs * (k2 * term1 * (ss * cfs - term1) + Gq * (k * css - q)) +
                    k2 * cf * s * ss * sfs * sfs)
        c22 = cs * (cfs * (k2 * term1 * (ss * cfs - term1) + Gqt * (k * css - q)) +
                    k2 * cf * s * ss * sfs * sfs)
        c31 = k * s * (term1 * cfs * (k * css - q) -
                       Gq * (cfs * (ss * cfs - term1) + ss * sfs * sfs))
        c32 = k * s * (term1 * cfs * (k * css - q) -
                       Gqt * (cfs * (ss * cfs - term1) - ss * sfs * sfs))
        c41 = k * cs * (cfs * css * (k * css - q) + k * ss * (ss * cfs - term1))
        c42 = k * cs * (cfs * css * (k * css - q) + k * ss * (ss * cfs - term1))
        c51 = Gq * (cfs * css * (q - k * css) - k * ss * (ss * cfs - term1))
        c52 = Gqt * (cfs * css * (q - k * css) - k * ss * (ss * cfs - term1))
    else:
        c21 = Gq * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                    k * s * ss * sfs * sfs)
        c22 = Gqt * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                     k * s * ss * sfs * sfs)
        c31 = k * ss * (k * cs * (ss * cfs - term1) + s * (kz + q))
        c32 = c31
        c41 = k * css * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                         k * s * ss * sfs * sfs)
        c42 = c41
        c51 = -css * (k2 * ss * (ss * cfs - term1) + Gq * cfs * (kz + q))
        c52 = -css * (k2 * ss * (ss * cfs - term1) + Gqt * cfs * (kz + q))

    qz = kz
    qzt = k * np.sqrt(er - s * s)

    vv = (1 + Rvi) * (-(1 - Rvi) * c11 / qz + (1 + Rvi) * c12 / qzt) + \
         (1 - Rvi) * ((1 - Rvi) * c21 / qz - (1 + Rvi) * c22 / qzt) + \
         (1 + Rvi) * ((1 - Rvi) * c31 / qz - (1 + Rvi) * c32 / (er * qzt)) + \
         (1 - Rvi) * ((1 + Rvi) * c41 / qz - er * (1 - Rvi) * c42 / qzt) + \
         (1 + Rvi) * ((1 + Rvi) * c51 / qz - (1 - Rvi) * c52 / qzt)

    hh = (1 + Rhi) * ((1 - Rhi) * c11 / qz - er * (1 + Rhi) * c12 / qzt) - \
         (1 - Rhi) * ((1 - Rhi) * c21 / qz - (1 + Rhi) * c22 / qzt) - \
         (1 + Rhi) * ((1 - Rhi) * c31 / qz - (1 + Rhi) * c32 / qzt) - \
         (1 - Rhi) * ((1 + Rhi) * c41 / qz - (1 - Rhi) * c42 / qzt) - \
         (1 + Rhi) * ((1 + Rhi) * c51 / qz - (1 - Rhi) * c52 / qzt)

    return vv, hh

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple I2EM_Bistat_model(double fr, double sig, double L, double thi,
                              double ths, double phs, complex er, int sp,
                              double xx, bint notify=False):
    if notify:
        print("I2EM_Bistat_model: using Cython backend")
    sig_cm = sig * 100.0
    L_cm = L * 100.0
    k = 2.0 * math.pi * fr / 30.0

    theta_i = math.radians(thi)
    theta_s = math.radians(ths)
    phi_s = math.radians(phs)
    phi_i = 0.0

    cs = math.cos(theta_i + 0.01)
    s = math.sin(theta_i + 0.01)
    ss = math.sin(theta_s)
    css = math.cos(theta_s)
    sf = math.sin(phi_i)
    cf = math.cos(phi_i)
    sfs = math.sin(phi_s)
    cfs = math.cos(phi_s)

    ks = k * sig_cm
    kl = k * L_cm
    ks2 = ks * ks

    kx = k * s * cf
    ky = k * s * sf
    kz = k * cs

    ksx = k * ss * cfs
    ksy = k * ss * sfs
    ksz = k * css

    s2 = s * s
    rt = np.sqrt(er - s2)
    Rvi = (er * cs - rt) / (er * cs + rt)
    Rhi = (cs - rt) / (cs + rt)

    wvnb = k * math.sqrt((ss * cfs - s * cf) ** 2 + (ss * sfs - s * sf) ** 2)

    Ts = 1
    error = 1e8
    while error > 1e-8:
        Ts += 1
        error = (ks2 * (cs + css) ** 2) ** Ts / factorial_double(Ts)

    wn, rss = roughness_spectrum(sp, xx, wvnb, sig_cm, L_cm, Ts)

    sqrt_er = np.sqrt(er)
    Rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
    Rh0 = -Rv0

    Ft = 8.0 * Rv0 * Rv0 * ss * (cs + np.sqrt(er - s2)) / (cs * np.sqrt(er - s2))
    a1 = 0.0
    b1 = 0.0
    for n in range(1, Ts + 1):
        a0 = (ks * cs) ** (2 * n) / factorial_double(n)
        wn_val = wn[n - 1]
        a1 += a0 * wn_val
        term = abs(Ft / 2.0 + 2.0 ** (n + 1) * Rv0 / cs * math.exp(-(ks * cs) ** 2)) ** 2
        b1 += a0 * term * wn_val

    St = 0.25 * abs(Ft) ** 2 * a1 / b1 if b1 != 0 else 0.0
    St0 = 1.0 / abs(1.0 + 8.0 * Rv0 / (cs * Ft)) ** 2
    Tf = 1.0 - St / St0 if St0 != 0 else 0.0

    sigx = 1.1 * sig_cm / L_cm
    sigy = sigx
    xxx = 3.0 * sigx

    Rav, Rah = _compute_average_reflectivities(cs, s, er, s2, sigx, sigy, xxx)

    if math.isclose(thi, ths) and math.isclose(phs, 180.0):
        Rvt = Rvi + (Rv0 - Rvi) * Tf
        Rht = Rhi + (Rh0 - Rhi) * Tf
    else:
        Rvt = Rav
        Rht = Rah

    fvv = 2.0 * Rvt * (s * ss - (1.0 + cs * css) * cfs) / (cs + css)
    fhh = -2.0 * Rht * (s * ss - (1.0 + cs * css) * cfs) / (cs + css)

    Fvvupi, Fhhupi = Fppupdn_is_calculations(+1, 1, Rvi, Rhi, er, k, kz, ksz,
                                             s, cs, ss, css, cf, cfs, sfs)
    Fvvups, Fhhups = Fppupdn_is_calculations(+1, 2, Rvi, Rhi, er, k, kz, ksz,
                                             s, cs, ss, css, cf, cfs, sfs)
    Fvvdni, Fhhdni = Fppupdn_is_calculations(-1, 1, Rvi, Rhi, er, k, kz, ksz,
                                             s, cs, ss, css, cf, cfs, sfs)
    Fvvdns, Fhhdns = Fppupdn_is_calculations(-1, 2, Rvi, Rhi, er, k, kz, ksz,
                                             s, cs, ss, css, cf, cfs, sfs)

    qi = k * cs
    qs = k * css

    Ivv = np.zeros(Ts, dtype=np.complex128)
    Ihh = np.zeros(Ts, dtype=np.complex128)

    for n in range(1, Ts + 1):
        idx = n - 1
        term_vv = 0.25 * (
            Fvvupi * (ksz - qi) ** (n - 1) * math.exp(-sig_cm ** 2 * (qi ** 2 - qi * (ksz - kz))) +
            Fvvdni * (ksz + qi) ** (n - 1) * math.exp(-sig_cm ** 2 * (qi ** 2 + qi * (ksz - kz))) +
            Fvvups * (kz + qs) ** (n - 1) * math.exp(-sig_cm ** 2 * (qs ** 2 - qs * (ksz - kz))) +
            Fvvdns * (kz - qs) ** (n - 1) * math.exp(-sig_cm ** 2 * (qs ** 2 + qs * (ksz - kz)))
        )
        Ivv[idx] = (kz + ksz) ** n * fvv * math.exp(-sig_cm ** 2 * kz * ksz) + term_vv

        term_hh = 0.25 * (
            Fhhupi * (ksz - qi) ** (n - 1) * math.exp(-sig_cm ** 2 * (qi ** 2 - qi * (ksz - kz))) +
            Fhhdni * (ksz + qi) ** (n - 1) * math.exp(-sig_cm ** 2 * (qi ** 2 + qi * (ksz - kz))) +
            Fhhups * (kz + qs) ** (n - 1) * math.exp(-sig_cm ** 2 * (qs ** 2 - qs * (ksz - kz))) +
            Fhhdns * (kz - qs) ** (n - 1) * math.exp(-sig_cm ** 2 * (qs ** 2 + qs * (ksz - kz)))
        )
        Ihh[idx] = (kz + ksz) ** n * fhh * math.exp(-sig_cm ** 2 * kz * ksz) + term_hh

    if math.isclose(thi, ths) and math.isclose(phs, 180.0):
        ct = 1.0 / math.tan(theta_i)
        cts = 1.0 / math.tan(theta_s)
        ctorslp = ct / math.sqrt(2.0) / rss
        ctsorslp = cts / math.sqrt(2.0) / rss
        shadf = 0.5 * (math.exp(-ctorslp ** 2) / (math.sqrt(math.pi) * ctorslp) - erfc(ctorslp))
        shadfs = 0.5 * (math.exp(-ctsorslp ** 2) / (math.sqrt(math.pi) * ctsorslp) - erfc(ctsorslp))
        ShdwS = 1.0 / (1.0 + shadf + shadfs)
    else:
        ShdwS = 1.0

    sigmavv = 0.0
    sigmahh = 0.0
    for n in range(1, Ts + 1):
        idx = n - 1
        a0 = wn[idx] / factorial_double(n) * sig_cm ** (2 * n)
        sigmavv += abs(Ivv[idx]) ** 2 * a0
        sigmahh += abs(Ihh[idx]) ** 2 * a0

    sigmavv *= ShdwS * k ** 2 / 2.0 * math.exp(-sig_cm ** 2 * (kz ** 2 + ksz ** 2))
    sigmahh *= ShdwS * k ** 2 / 2.0 * math.exp(-sig_cm ** 2 * (kz ** 2 + ksz ** 2))

    sigma_0_vv = 10.0 * math.log10(sigmavv)
    sigma_0_hh = 10.0 * math.log10(sigmahh)

    depol_ratio = 0.05
    sigma_0_hv = 10.0 * math.log10(depol_ratio * math.sqrt(sigmavv * sigmahh))
    sigma_0_vh = sigma_0_hv

    return sigma_0_vv, sigma_0_hh, sigma_0_hv, sigma_0_vh
