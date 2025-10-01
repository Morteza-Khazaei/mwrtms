import os
import math
import numpy as np
from functools import partial
from multiprocessing import get_context
from typing import Iterable, Sequence, Tuple, Optional

from scipy.special import erfc, gamma, kv as besselk, jv as besselj, factorial
from scipy.integrate import dblquad, quad

try:
    from ._i2em_cy import I2EM_Bistat_model as _I2EM_Bistat_model_cython
except ImportError:  # pragma: no cover - fallback only exercised without extension
    _I2EM_Bistat_model_cython = None

def _I2EM_Bistat_model_python(fr, sig, L, thi, ths, phs, er, sp, xx, notify=False):
    """
    I2EM bistatic scattering model for single-scale random surfaces.

    Parameters:
    - fr   : frequency (GHz)
    - sig  : RMS height (m)
    - L    : correlation length (m)
    - thi  : incidence angle (deg)
    - ths  : scattering angle (deg)
    - phs  : relative azimuth angle (deg)
    - er   : complex dielectric constant
    - sp   : correlation function type (1-exp, 2-gauss, 3-power, 4-x-exp)
    - xx   : parameter for sp=3,4

    Returns:
    - sigma_0_vv, sigma_0_hh, sigma_0_hv, sigma_0_vh (dB)

    Examples
    --------
    >>> I2EM_Bistat_model(
    ...     fr=10.0,
    ...     sig=0.01,
    ...     L=0.1,
    ...     thi=30.0,
    ...     ths=30.0,
    ...     phs=180.0,
    ...     er=3.0 + 0.5j,
    ...     sp=1,
    ...     xx=0.0,
    ...     notify=False,
    ... )  # doctest: +SKIP
    """

    # --- Helper functions ---

    def roughness_spectrum(sp, xx, wvnb, sig, L, Ts):
        wn = np.zeros(Ts)
        if sp == 1:  # exponential
            for n in range(1, Ts + 1):
                wn[n - 1] = (L**2) / n**2 * (1 + (wvnb * L / n)**2)**(-1.5)
            rss = sig / L
        elif sp == 2:  # gaussian
            for n in range(1, Ts + 1):
                wn[n - 1] = L**2 / (2 * n) * np.exp(-(wvnb * L)**2 / (4 * n))
            rss = np.sqrt(2) * sig / L
        elif sp == 3:  # x-power
            for n in range(1, Ts + 1):
                if wvnb == 0:
                    wn[n - 1] = L**2 / (3 * n - 2)
                else:
                    arg = wvnb * L
                    nu = 1 - xx * n
                    wn[n - 1] = (
                        L**2 * arg**(-nu) * besselk(nu, arg) /
                        (2**(xx * n - 1) * gamma(xx * n))
                    )
            rss = np.sqrt(2 * xx) * sig / L
        elif sp == 4:  # x-exponential
            def integrand(z, L=L, wvnb=wvnb, n=1, xx=xx):
                return np.exp(-np.abs(z)**xx) * besselj(0, z * wvnb * L / n**(1 / xx)) * z
            for n in range(1, Ts + 1):
                integral, _ = quad(integrand, 0, 9)
                wn[n - 1] = L**2 / n**(2 / xx) * integral
            rss = sig / L
        else:
            raise ValueError("Unsupported sp value, must be 1-4")
        return wn, rss

    def Rav_integration(Zx, Zy, cs, s, er, s2, sigx, sigy):
        A = cs + Zx * s
        B = er * (1 + Zx**2 + Zy**2)
        CC = s2 - 2 * Zx * s * cs + Zx**2 * cs**2 + Zy**2
        Rv = (er * A - np.sqrt(B - CC)) / (er * A + np.sqrt(B - CC))
        pd = np.exp(-Zx**2 / (2 * sigx**2) - Zy**2 / (2 * sigy**2))
        return Rv * pd

    def Rah_integration(Zx, Zy, cs, s, er, s2, sigx, sigy):
        A = cs + Zx * s
        B = er * (1 + Zx**2 + Zy**2)
        CC = s2 - 2 * Zx * s * cs + Zx**2 * cs**2 + Zy**2
        Rh = (A - np.sqrt(B - CC)) / (A + np.sqrt(B - CC))
        pd = np.exp(-Zx**2 / (2 * sigx**2) - Zy**2 / (2 * sigy**2))
        return Rh * pd

    def Fppupdn_is_calculations(ud, is_mode, Rvi, Rhi, er, k, kz, ksz,
                                s, cs, ss, css, cf, cfs, sfs):
        if is_mode == 1:
            Gq = ud * kz
            Gqt = ud * k * np.sqrt(er - s**2)
            q = ud * kz
        elif is_mode == 2:
            Gq = ud * ksz
            Gqt = ud * k * np.sqrt(er - ss**2)
            q = ud * ksz
        else:
            raise ValueError("is_mode must be 1 or 2")

        k2 = k**2
        term1 = s * cf
        term2 = ss * cfs
        term3 = ss * sfs

        c11 = k * cfs * (ksz - q)
        c12 = k * cfs * (ksz - q)

        if is_mode == 1:
            c21 = cs * (cfs * (k2 * term1 * (ss * cfs - term1) + Gq * (k * css - q)) +
                        k2 * cf * s * ss * sfs**2)
            c22 = cs * (cfs * (k2 * term1 * (ss * cfs - term1) + Gqt * (k * css - q)) +
                        k2 * cf * s * ss * sfs**2)
            c31 = k * s * (term1 * cfs * (k * css - q) -
                           Gq * (cfs * (ss * cfs - term1) + ss * sfs**2))
            c32 = k * s * (term1 * cfs * (k * css - q) -
                           Gqt * (cfs * (ss * cfs - term1) - ss * sfs**2))
            c41 = k * cs * (cfs * css * (k * css - q) + k * ss * (ss * cfs - term1))
            c42 = k * cs * (cfs * css * (k * css - q) + k * ss * (ss * cfs - term1))
            c51 = Gq * (cfs * css * (q - k * css) - k * ss * (ss * cfs - term1))
            c52 = Gqt * (cfs * css * (q - k * css) - k * ss * (ss * cfs - term1))
        else:
            c21 = Gq * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                        k * s * ss * sfs**2)
            c22 = Gqt * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                         k * s * ss * sfs**2)
            c31 = k * ss * (k * cs * (ss * cfs - term1) + s * (kz + q))
            c32 = c31
            c41 = k * css * (cfs * (cs * (kz + q) - k * s * (ss * cfs - term1)) -
                             k * s * ss * sfs**2)
            c42 = c41
            c51 = -css * (k2 * ss * (ss * cfs - term1) + Gq * cfs * (kz + q))
            c52 = -css * (k2 * ss * (ss * cfs - term1) + Gqt * cfs * (kz + q))

        qz = kz
        qzt = k * np.sqrt(er - s**2)

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

    # --- Begin main calculation ---

    # Convert units: meters to cm
    sig_cm = sig * 100
    L_cm = L * 100
    k = 2 * np.pi * fr / 30  # wavenumber (cm^-1)

    # Angles in radians
    theta_i = np.radians(thi)
    theta_s = np.radians(ths)
    phi_s = np.radians(phs)
    phi_i = 0.0

    # Trig terms
    cs = np.cos(theta_i + 0.01)
    s = np.sin(theta_i + 0.01)
    ss = np.sin(theta_s)
    css = np.cos(theta_s)
    sf = np.sin(phi_i)
    cf = np.cos(phi_i)
    sfs = np.sin(phi_s)
    cfs = np.cos(phi_s)

    ks = k * sig_cm
    kl = k * L_cm
    ks2 = ks ** 2

    kx = k * s * cf
    ky = k * s * sf
    kz = k * cs

    ksx = k * ss * cfs
    ksy = k * ss * sfs
    ksz = k * css

    s2 = s**2
    rt = np.sqrt(er - s2)
    Rvi = (er * cs - rt) / (er * cs + rt)
    Rhi = (cs - rt) / (cs + rt)

    wvnb = k * np.sqrt((ss * cfs - s * cf)**2 + (ss * sfs - s * sf)**2)

    # Series order convergence
    Ts = 1
    error = 1e8
    while error > 1e-8:
        Ts += 1
        error = (ks2 * (cs + css)**2) ** Ts / math.factorial(Ts)

    # Surface roughness spectrum
    wn, rss = roughness_spectrum(sp, xx, wvnb, sig_cm, L_cm, Ts)

    # Fresnel transition terms
    Rv0 = (np.sqrt(er) - 1) / (np.sqrt(er) + 1)
    Rh0 = -Rv0

    Ft = 8 * Rv0**2 * ss * (cs + np.sqrt(er - s2)) / (cs * np.sqrt(er - s2))
    a1 = 0
    b1 = 0
    for n in range(1, Ts + 1):
        a0 = (ks * cs) ** (2 * n) / math.factorial(n)
        a1 += a0 * wn[n - 1]
        term = np.abs(Ft / 2 + 2 ** (n + 1) * Rv0 / cs * np.exp(-(ks * cs) ** 2)) ** 2
        b1 += a0 * term * wn[n - 1]

    St = 0.25 * np.abs(Ft) ** 2 * a1 / b1
    St0 = 1 / np.abs(1 + 8 * Rv0 / (cs * Ft)) ** 2
    Tf = 1 - St / St0

    # Average reflectivities using slope PDF
    sigx = 1.1 * sig_cm / L_cm
    sigy = sigx
    xxx = 3 * sigx

    Rav, _ = dblquad(
        lambda Zy, Zx: Rav_integration(Zx, Zy, cs, s, er, s2, sigx, sigy),
        -xxx, xxx, lambda x: -xxx, lambda x: xxx
    )
    Rah, _ = dblquad(
        lambda Zy, Zx: Rah_integration(Zx, Zy, cs, s, er, s2, sigx, sigy),
        -xxx, xxx, lambda x: -xxx, lambda x: xxx
    )

    Rav /= 2 * np.pi * sigx * sigy
    Rah /= 2 * np.pi * sigx * sigy

    # Select reflection coefficients
    if np.isclose(thi, ths) and np.isclose(phs, 180):
        Rvt = Rvi + (Rv0 - Rvi) * Tf
        Rht = Rhi + (Rh0 - Rhi) * Tf
    else:
        Rvt = Rav
        Rht = Rah

    # Scattering amplitude terms
    fvv = 2 * Rvt * (s * ss - (1 + cs * css) * cfs) / (cs + css)
    fhh = -2 * Rht * (s * ss - (1 + cs * css) * cfs) / (cs + css)

    # Calculate Fppupdn terms
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

    Ivv = np.zeros(Ts)
    Ihh = np.zeros(Ts)

    for n in range(1, Ts + 1):
        idx = n - 1
        term_vv = (
            0.25 * (
                Fvvupi * (ksz - qi) ** (n - 1) * np.exp(-sig_cm**2 * (qi**2 - qi * (ksz - kz))) +
                Fvvdni * (ksz + qi) ** (n - 1) * np.exp(-sig_cm**2 * (qi**2 + qi * (ksz - kz))) +
                Fvvups * (kz + qs) ** (n - 1) * np.exp(-sig_cm**2 * (qs**2 - qs * (ksz - kz))) +
                Fvvdns * (kz - qs) ** (n - 1) * np.exp(-sig_cm**2 * (qs**2 + qs * (ksz - kz)))
            )
        )
        Ivv[idx] = (kz + ksz) ** n * fvv * np.exp(-sig_cm**2 * kz * ksz) + term_vv

        term_hh = (
            0.25 * (
                Fhhupi * (ksz - qi) ** (n - 1) * np.exp(-sig_cm**2 * (qi**2 - qi * (ksz - kz))) +
                Fhhdni * (ksz + qi) ** (n - 1) * np.exp(-sig_cm**2 * (qi**2 + qi * (ksz - kz))) +
                Fhhups * (kz + qs) ** (n - 1) * np.exp(-sig_cm**2 * (qs**2 - qs * (ksz - kz))) +
                Fhhdns * (kz - qs) ** (n - 1) * np.exp(-sig_cm**2 * (qs**2 + qs * (ksz - kz)))
            )
        )
        Ihh[idx] = (kz + ksz) ** n * fhh * np.exp(-sig_cm**2 * kz * ksz) + term_hh

    # Shadowing function
    if np.isclose(thi, ths) and np.isclose(phs, 180):
        ct = 1 / np.tan(theta_i)
        cts = 1 / np.tan(theta_s)
        ctorslp = ct / np.sqrt(2) / rss
        ctsorslp = cts / np.sqrt(2) / rss
        shadf = 0.5 * (np.exp(-ctorslp**2) / (np.sqrt(np.pi) * ctorslp) - erfc(ctorslp))
        shadfs = 0.5 * (np.exp(-ctsorslp**2) / (np.sqrt(np.pi) * ctsorslp) - erfc(ctsorslp))
        ShdwS = 1 / (1 + shadf + shadfs)
    else:
        ShdwS = 1

    # Final sigma calculations
    sigmavv = 0
    sigmahh = 0
    for n in range(1, Ts + 1):
        idx = n - 1
        a0 = wn[idx] / math.factorial(n) * sig_cm ** (2 * n)
        sigmavv += np.abs(Ivv[idx]) ** 2 * a0
        sigmahh += np.abs(Ihh[idx]) ** 2 * a0

    sigmavv *= ShdwS * k**2 / 2 * np.exp(-sig_cm**2 * (kz**2 + ksz**2))
    sigmahh *= ShdwS * k**2 / 2 * np.exp(-sig_cm**2 * (kz**2 + ksz**2))

    sigma_0_vv = 10 * np.log10(sigmavv)
    sigma_0_hh = 10 * np.log10(sigmahh)

    # Empirical cross-pol approximation
    depol_ratio = 0.05
    sigma_0_hv = 10 * np.log10(depol_ratio * np.sqrt(sigmavv * sigmahh))
    sigma_0_vh = sigma_0_hv

    if notify:
        print("I2EM_Bistat_model: using Python backend")

    return sigma_0_vv, sigma_0_hh, sigma_0_hv, sigma_0_vh


if _I2EM_Bistat_model_cython is not None:
    I2EM_Bistat_model = _I2EM_Bistat_model_cython
    I2EM_BACKEND = "cython"
else:
    I2EM_Bistat_model = _I2EM_Bistat_model_python
    I2EM_BACKEND = "python"


def _resolve_backend(backend: str):
    backend = backend.lower()
    if backend == "auto":
        func = I2EM_Bistat_model
        name = I2EM_BACKEND
    elif backend == "cython":
        if _I2EM_Bistat_model_cython is None:
            raise RuntimeError("Cython backend is unavailable. Rebuild the extension or install build dependencies.")
        func = _I2EM_Bistat_model_cython
        name = "cython"
    elif backend == "python":
        func = _I2EM_Bistat_model_python
        name = "python"
    else:
        raise ValueError("backend must be one of {'auto', 'cython', 'python'}")
    return func, name


def _i2em_call(args, backend: str):
    func, _ = _resolve_backend(backend)
    return func(*args, notify=False)


def get_i2em_backend() -> str:
    """Return the active backend identifier."""
    return I2EM_BACKEND


def I2EM_Bistat_model_batch(
    parameters: Iterable[Sequence],
    *,
    backend: str = "auto",
    processes: Optional[int] = os.cpu_count()/3,
    chunksize: int = 1,
    mp_start_method: str = "spawn",
    notify_backend: bool = False,
) -> Tuple[Tuple[float, float, float, float], ...]:
    """Evaluate I2EM for many parameter tuples in parallel.

    Parameters
    ----------
    parameters:
        Iterable of argument sequences. Each element must match the positional
        signature ``(fr, sig, L, thi, ths, phs, er, sp, xx)``.
    backend:
        'auto' (default) chooses the active backend, 'cython' forces the
        compiled version, and 'python' forces the pure-Python fallback.
    processes:
        Number of worker processes to use. Defaults to ``os.cpu_count()``.
    chunksize:
        Work chunk size forwarded to ``ProcessPoolExecutor.map``.
    mp_start_method:
        Multiprocessing start method. Use 'spawn' by default for safety with
        C extensions; set to 'fork' on POSIX only if you understand the
        implications.
    notify_backend:
        When ``True`` prints the backend that will be used before dispatching.

    Returns
    -------
    tuple
        Tuple of I2EM results in the same order as ``parameters``.
    """

    _, backend_name = _resolve_backend(backend)

    # Eagerly materialise parameters for repeated iteration and validation
    params_list = [tuple(p) for p in parameters]
    if not params_list:
        return tuple()

    for idx, item in enumerate(params_list):
        if len(item) != 9:
            raise ValueError(
                f"Parameter tuple at index {idx} must contain 9 elements (fr, sig, L, thi, ths, phs, er, sp, xx); got {len(item)}"
            )

    process_label = processes if processes is not None else "default"
    if notify_backend:
        print(f"I2EM_Bistat_model_batch: using {backend_name} backend across {process_label} processes")

    if processes == 1:
        results = [_i2em_call(args, backend) for args in params_list]
    else:
        worker = partial(_i2em_call, backend=backend)
        with get_context(mp_start_method).Pool(processes=processes) as pool:
            results = pool.map(worker, params_list, chunksize)

    return tuple(results)


__all__ = [
    "I2EM_Bistat_model",
    "_I2EM_Bistat_model_python",
    "I2EM_BACKEND",
    "get_i2em_backend",
    "I2EM_Bistat_model_batch",
]
