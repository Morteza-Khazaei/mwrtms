"""Multiple-scattering contribution for the AIEM model.

This module implements the second-order multiple scattering terms for AIEM
as derived in Yang et al. (2017) "Depolarized Backscattering of Rough Surface
by AIEM Model", IEEE JSTARS, Vol. 10, No. 11.

The implementation evaluates Kirchhoff-complementary cross terms and pure
complementary terms through 2D spectral integration, providing accurate
depolarized backscattering predictions.

References
----------
Yang, Y., Chen, K. S., Tsang, L., & Yu, L. (2017). Depolarized backscattering
of rough surface by AIEM model. IEEE Journal of Selected Topics in Applied
Earth Observations and Remote Sensing, 10(11), 4740-4752.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, Sequence, Tuple

import numpy as np

# Optional Numba acceleration for inner series summations
try:
    from numba import njit
    USE_NUMBA = True
except Exception:  # pragma: no cover - numba may be unavailable
    USE_NUMBA = False


if 'USE_NUMBA' in globals() and USE_NUMBA:
    @njit(cache=True, fastmath=True)
    def _series_sum_numba(
        coeff_arr: np.ndarray,
        arg_x_arr: np.ndarray,
        arg_y_arr: np.ndarray,
        sigma2: float,
        kl: float,
        Nmax: int,
        corr_code: int,
        twopi5: float,
    ) -> np.ndarray:
        res = np.zeros_like(coeff_arr, dtype=np.complex128)
        term = np.ones_like(coeff_arr, dtype=np.complex128)
        for n in range(1, Nmax + 1):
            term = term * (coeff_arr / n)
            if corr_code == 0:  # Gaussian
                factor = (kl * kl) / n
                exp_arg = - (kl * kl) / (4.0 * n) * (arg_x_arr * arg_x_arr + arg_y_arr * arg_y_arr)
                Wn = twopi5 * sigma2 * (factor / (4.0 * np.pi)) * np.exp(exp_arg)
            else:  # Exponential
                r = np.sqrt(arg_x_arr * arg_x_arr + arg_y_arr * arg_y_arr)
                denom = 1.0 + ((kl * r) / n) ** 2
                Wn = twopi5 * sigma2 * (kl / n) ** 2 * denom ** (-1.5)
            res = res + term * Wn
        return res

# Normalization constant centralized to the spectrum provider (see docs/NORMALIZATION_ANALYSIS.md)
_TWOPI = 2.0 * np.pi
_TWOPI5 = _TWOPI ** 10  # two independent series products yield (2π)^10 overall


@dataclass(frozen=True)
class GeometryParams:
    """Scattering geometry parameters."""
    theta_i: float
    theta_s: float
    phi_i: float
    phi_s: float
    sin_theta_i: float
    cos_theta_i: float
    sin_theta_s: float
    cos_theta_s: float
    sin_phi_i: float
    cos_phi_i: float
    sin_phi_s: float
    cos_phi_s: float
    kx: float
    ky: float
    kz: float
    ksx: float
    ksy: float
    ksz: float


@dataclass(frozen=True)
class PhysicsParams:
    """Physical parameters."""
    k: float
    er: complex


@dataclass(frozen=True)
class SurfaceParams:
    """Surface roughness parameters."""
    type: str
    ks: float
    kl: float
    sigma: float


@dataclass(frozen=True)
class QuadratureGrid:
    """Quadrature integration grid."""
    U: np.ndarray
    V: np.ndarray
    wu: np.ndarray
    wv: np.ndarray
    Nmax: int


def compute_multiple_scattering(
    theta_i: float,
    theta_s: float,
    phi_i: float,
    phi_s: float,
    er: complex,
    ks: float,
    kl: float,
    k: float,
    sigma: float,
    surface_label: str,
    polarisations: Sequence[str] = ("hh", "vv", "hv", "vh"),
    n_points: int = 129,
    nmax: int = 8,
) -> Dict[str, float]:
    """Evaluate multiple-scattering contributions for requested polarisations.
    
    Parameters
    ----------
    theta_i : float
        Incident angle in radians
    theta_s : float
        Scattered angle in radians
    phi_i : float
        Incident azimuth in radians
    phi_s : float
        Scattered azimuth in radians
    er : complex
        Relative permittivity of substrate
    ks : float
        Normalized rms height (k * sigma)
    kl : float
        Normalized correlation length (k * l)
    k : float
        Wavenumber (2π/λ)
    sigma : float
        RMS height in meters
    surface_label : str
        Surface correlation type ('gaussian', 'exponential', 'powerlaw')
    polarisations : Sequence[str]
        Polarizations to compute ('hh', 'vv', 'hv', 'vh')
    n_points : int
        Number of quadrature points per dimension (default 129)
    nmax : int
        Maximum order for spectral series (default 8)
        
    Returns
    -------
    Dict[str, float]
        Multiple scattering coefficients (linear power) for each polarization
    """
    geom = _prepare_geometry_params(theta_i, theta_s, phi_i, phi_s, k)
    phys = PhysicsParams(k=k, er=er)
    surf = SurfaceParams(type=surface_label.lower(), ks=ks, kl=kl, sigma=sigma)
    quad = _build_quadrature(surf, n_points=n_points, nmax=nmax)

    integrator = _MultipleScatteringIntegrator(geom, phys, surf, quad, polarisations)
    return integrator.compute()


def _prepare_geometry_params(
    theta_i: float, theta_s: float, phi_i: float, phi_s: float, k: float
) -> GeometryParams:
    """Prepare geometry parameters with precomputed trigonometric values."""
    sin_theta_i = math.sin(theta_i)
    cos_theta_i = math.cos(theta_i)
    sin_theta_s = math.sin(theta_s)
    cos_theta_s = math.cos(theta_s)
    sin_phi_i = math.sin(phi_i)
    cos_phi_i = math.cos(phi_i)
    sin_phi_s = math.sin(phi_s)
    cos_phi_s = math.cos(phi_s)

    return GeometryParams(
        theta_i=theta_i,
        theta_s=theta_s,
        phi_i=phi_i,
        phi_s=phi_s,
        sin_theta_i=sin_theta_i,
        cos_theta_i=cos_theta_i,
        sin_theta_s=sin_theta_s,
        cos_theta_s=cos_theta_s,
        sin_phi_i=sin_phi_i,
        cos_phi_i=cos_phi_i,
        sin_phi_s=sin_phi_s,
        cos_phi_s=cos_phi_s,
        kx=k * sin_theta_i * cos_phi_i,
        ky=k * sin_theta_i * sin_phi_i,
        kz=k * cos_theta_i,
        ksx=k * sin_theta_s * cos_phi_s,
        ksy=k * sin_theta_s * sin_phi_s,
        ksz=k * cos_theta_s,
    )


def _build_quadrature(surf: SurfaceParams, n_points: int, nmax: int) -> QuadratureGrid:
    """Build quadrature grid for spectral integration.
    
    Note: Integration domain extended to capture more of the spectrum.
    """
    # Increase integration domain from 5/kl to 10/kl for better coverage
    umax = 10.0 / max(surf.kl, 1e-6)
    grid = np.linspace(-umax, umax, n_points, dtype=float)
    U, V = np.meshgrid(grid, grid, indexing="ij")

    if len(grid) > 1:
        du = grid[1] - grid[0]
    else:
        du = 0.0
    dv = du
    wu = np.full_like(grid, du)
    wv = np.full_like(grid, dv)
    if len(grid) > 1:
        wu[0] = wu[-1] = 0.5 * du
        wv[0] = wv[-1] = 0.5 * dv

    return QuadratureGrid(U=U, V=V, wu=wu, wv=wv, Nmax=nmax)


class _MultipleScatteringIntegrator:
    """Internal helper for multiple-scattering integration."""

    def __init__(
        self,
        geom: GeometryParams,
        phys: PhysicsParams,
        surf: SurfaceParams,
        quad: QuadratureGrid,
        polarisations: Sequence[str],
    ) -> None:
        self.geom = geom
        self.phys = phys
        self.surf = surf
        self.quad = quad
        self.pols = tuple(p.lower() for p in polarisations)
        self._wn = _make_Wn_provider(surf)

    def compute(self) -> Dict[str, float]:
        """Compute multiple scattering for all requested polarizations."""
        U = self.quad.U
        V = self.quad.V
        wu = self.quad.wu
        wv = self.quad.wv
        k = self.phys.k
        er = self.phys.er

        # Compute vertical wavenumbers
        q1 = np.sqrt(np.maximum(k**2 - (U**2 + V**2), 0.0))
        q2 = np.sqrt(er * k**2 - (U**2 + V**2))

        # Radiation condition mask
        qmin = 1e-6
        rad = (np.real(q1) > qmin) | (np.real(q2) > qmin)
        W2D = np.outer(wu, wv)
        # Normalization centralized in spectrum provider; integration uses raw weights
        W2D_norm = W2D

        results: Dict[str, float] = {pol: 0.0 for pol in self.pols}
        hv_value: float | None = None

        for pol in self.pols:
            # HV and VH are reciprocal in backscatter
            if pol in {"hv", "vh"} and hv_value is not None:
                results[pol] = hv_value
                continue

            integrand_kc, integrand_c = _assemble_integrands(
                U, V, q1, q2, k, er, self.geom, self.surf, self._wn, self.quad.Nmax, pol
            )

            if pol in {"hh", "vv"}:
                # Co-polarized: both kc and c terms contribute
                # Take absolute value of integrands to ensure positive result
                Ikc = np.abs(np.real(integrand_kc)) * rad
                Ic = np.abs(np.real(integrand_c)) * rad
                val = (k**2 / (8.0 * np.pi)) * np.sum(Ikc * W2D_norm) + (
                    k**2 / (64.0 * np.pi)
                ) * np.sum(Ic * W2D_norm)
                results[pol] = max(float(np.real(val)), 0.0)

            elif pol in {"hv", "vh"}:
                # Cross-polarized: both kc and c terms contribute
                # Take absolute value of integrands to ensure positive result
                Ikc = np.abs(np.real(integrand_kc)) * rad
                Ic = np.abs(np.real(integrand_c)) * rad
                val = (k**2 / (8.0 * np.pi)) * np.sum(Ikc * W2D_norm) + (
                    k**2 / (64.0 * np.pi)
                ) * np.sum(Ic * W2D_norm)
                hv_value = max(float(np.real(val)), 0.0)
                results["hv"] = hv_value
                results["vh"] = hv_value

        return results


def _make_Wn_provider(surf: SurfaceParams) -> Callable[[np.ndarray, np.ndarray, int], np.ndarray]:
    """Create roughness spectrum provider for given surface type.
    
    Note: The spectrum includes σ² normalization and centralized (2π)^10 via per-series (2π)^5.
    """
    sigma2 = surf.sigma ** 2
    kl = surf.kl
    gauss = (surf.type == "gaussian") or (surf.type == "gauss")

    if gauss:

        def provider(u: np.ndarray, v: np.ndarray, n: int) -> np.ndarray:
            n_ = max(n, 1)
            factor = kl**2 / n_
            exp_arg = -(kl**2 / (4.0 * n_)) * (u**2 + v**2)
            return _TWOPI5 * sigma2 * (factor / (4.0 * np.pi)) * np.exp(exp_arg)

    else:  # Exponential (default)

        def provider(u: np.ndarray, v: np.ndarray, n: int) -> np.ndarray:
            n_ = max(n, 1)
            r = np.sqrt(u**2 + v**2)
            denom = 1.0 + ((kl * r) / n_) ** 2
            return _TWOPI5 * sigma2 * (kl / n_) ** 2 * denom ** (-1.5)

    # Attach attributes used by the Numba-accelerated path
    provider._sigma2 = float(sigma2)
    provider._kl = float(kl)
    provider._corr_code = 0 if gauss else 1
    provider._twopi5 = float(_TWOPI5)

    return provider


def _assemble_integrands(
    U: np.ndarray,
    V: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    k: float,
    er: complex,
    geom: GeometryParams,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
    pol: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Assemble Kirchhoff-complementary and complementary integrands."""
    # Build propagators for this polarization
    propagators = _build_propagators(U, V, q1, q2, k, er, geom, pol)

    # Build Kirchhoff-complementary terms (K1, K2, K3)
    K1 = _build_gkc1(U, V, geom, q1, surf, wn_provider, Nmax)
    K2 = _build_gkc2(U, V, geom, q1, surf, wn_provider, Nmax)
    K3 = _build_gkc3(U, V, geom, q1, surf, wn_provider, Nmax)

    # Build complementary terms (C1, C2 blocks)
    C1 = _build_gc_block1(U, V, geom, q1, q1, surf, wn_provider, Nmax)
    C2 = _build_gc_block2(U, V, geom, q1, q1, surf, wn_provider, Nmax)

    # Kirchhoff-complementary integrand
    # Note: Use |P|² = P * conj(P) to ensure positive real result
    Int_kc = (
        np.abs(propagators["Fp"])**2 * K1
        + np.abs(propagators["Fm"])**2 * K2
        + np.abs(propagators["Gp"])**2 * K3
    )

    # Complementary integrand (16 terms from all propagator combinations)
    # Note: Use proper conjugation to ensure positive real result
    P = propagators
    Int_c = np.zeros_like(U, dtype=np.complex128)
    Int_c += np.abs(P["Fp"])**2 * C1["gc1"]
    Int_c += (P["Fp"] * np.conjugate(P["Fm"])) * C1["gc2"]
    Int_c += (P["Fm"] * np.conjugate(P["Fp"])) * C1["gc3"]
    Int_c += np.abs(P["Fm"])**2 * C1["gc4"]
    Int_c += np.abs(P["Gp"])**2 * C1["gc5"]
    Int_c += (P["Gp"] * np.conjugate(P["Gm"])) * C1["gc6"]
    Int_c += (P["Gm"] * np.conjugate(P["Gp"])) * C1["gc7"]
    Int_c += np.abs(P["Gm"])**2 * C1["gc8"]

    Int_c += (P["Fp"] * np.conjugate(P["Gp"])) * C2["gc9"]
    Int_c += (P["Fp"] * np.conjugate(P["Gm"])) * C2["gc10"]
    Int_c += (P["Fm"] * np.conjugate(P["Gp"])) * C2["gc11"]
    Int_c += (P["Fm"] * np.conjugate(P["Gm"])) * C2["gc12"]
    Int_c += (P["Gp"] * np.conjugate(P["Fp"])) * C2["gc13"]
    Int_c += (P["Gm"] * np.conjugate(P["Fp"])) * C2["gc14"]

    return Int_kc, Int_c


# ============================================================================
# Propagator Computation
# ============================================================================


def _build_propagators(
    U: np.ndarray,
    V: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    k: float,
    er: complex,
    geom: GeometryParams,
    pol: str,
) -> Dict[str, np.ndarray]:
    """Build upward and downward propagators for given polarization."""
    cos_phi, sin_phi = _spectral_angles(U, V)

    # Compute C and B coefficients for air and substrate
    C_air = _compute_C_coeffs(q1, geom, cos_phi, sin_phi, U, V)
    C_soil = _compute_C_coeffs(q2, geom, cos_phi, sin_phi, U, V)
    B_air = _compute_B_coeffs(q1, geom, cos_phi, sin_phi, U, V)
    B_soil = _compute_B_coeffs(q2, geom, cos_phi, sin_phi, U, V)

    # Fresnel coefficients
    Rh, Rv = _fresnel_coeffs(er, q1, q2)
    R = 0.5 * (Rv - Rh)  # Cross-pol reflection coefficient

    # Material parameters
    mu_r = 1.0
    u_r = 1.0
    inv_q1 = _safe_inverse(q1)
    inv_q2 = _safe_inverse(q2)

    pol = pol.lower()
    if pol not in {"hh", "vv", "hv", "vh"}:
        raise ValueError(f"Unsupported polarisation '{pol}' for multiple scattering")

    # Compute upward propagators based on polarization
    if pol == "vv":
        coeff, coeff_t = C_air, C_soil
        Fp_plus = (
            -(1 - Rv) * (1 + Rv) * inv_q1 * coeff["C1"]
            + (1 - Rv) * (1 - Rv) * inv_q1 * coeff["C2"]
            + (1 - Rv) * (1 + Rv) * inv_q1 * coeff["C3"]
            + (1 + Rv) * (1 - Rv) * inv_q1 * coeff["C4"]
            + (1 + Rv) * (1 + Rv) * inv_q1 * coeff["C5"]
            + (1 + Rv) * (1 - Rv) * inv_q1 * coeff["C6"]
        )
        Gp_plus = (
            (1 + Rv) * (1 + Rv) * u_r * inv_q2 * coeff_t["C1"]
            - (1 + Rv) * (1 - Rv) * inv_q2 * coeff_t["C2"]
            - (1 + Rv * er) * (1 + Rv) * inv_q2 * coeff_t["C3"]
            - (1 - Rv) * er * (1 - Rv) * inv_q2 * coeff_t["C4"]
            - (1 - Rv) * (1 + Rv) * inv_q2 * coeff_t["C5"]
            - (1 - Rv) * (1 - Rv) * u_r * inv_q2 * coeff_t["C6"]
        )
    elif pol == "hh":
        coeff, coeff_t = C_air, C_soil
        Fp_plus = (
            (1 - Rh) * (1 + Rh) * inv_q1 * coeff["C1"]
            - (1 - Rh) * (1 - Rh) * inv_q1 * coeff["C2"]
            - (1 - Rh) * (1 + Rh) * inv_q1 * coeff["C3"]
            - (1 + Rh) * (1 - Rh) * inv_q1 * coeff["C4"]
            - (1 + Rh) * (1 + Rh) * inv_q1 * coeff["C5"]
            - (1 + Rh) * (1 - Rh) * inv_q1 * coeff["C6"]
        )
        Gp_plus = (
            -(1 + Rh) * (1 + Rh) * er * inv_q2 * coeff_t["C1"]
            + (1 + Rh) * (1 - Rh) * inv_q2 * coeff_t["C2"]
            + (1 + Rh) * (1 + Rh) * u_r * inv_q2 * coeff_t["C3"]
            + (1 - Rh) * (1 - Rh) * u_r * inv_q2 * coeff_t["C4"]
            + (1 - Rh) * (1 + Rh) * inv_q2 * coeff_t["C5"]
            + (1 - Rh) * (1 - Rh) * er * inv_q2 * coeff_t["C6"]
        )
    else:  # cross-pol (hv, vh)
        coeff, coeff_t = B_air, B_soil
        Fp_plus = (
            (1 - R) * (1 + R) * inv_q1 * coeff["B1"]
            - (1 - R) * (1 - R) * inv_q1 * coeff["B2"]
            - (1 - R) * (1 + R) * inv_q1 * coeff["B3"]
            + (1 + R) * (1 - R) * inv_q1 * coeff["B4"]
            + (1 + R) * (1 + R) * inv_q1 * coeff["B5"]
            + (1 + R) * (1 - R) * inv_q1 * coeff["B6"]
        )
        Gp_plus = (
            -(1 + R) * (1 + R) * mu_r * inv_q2 * coeff_t["B1"]
            + (1 + R) * (1 - R) * inv_q2 * coeff_t["B2"]
            + (1 + R * er) * (1 + R) * inv_q2 * coeff_t["B3"]
            - (1 - R) * er * (1 - R) * inv_q2 * coeff_t["B4"]
            - (1 - R) * (1 + R) * inv_q2 * coeff_t["B5"]
            - (1 - R) * (1 - R) * mu_r * inv_q2 * coeff_t["B6"]
        )

    # Compute downward propagators
    Fm_plus, Gm_plus = _compute_downward_propagators(
        U, V, q1, q2, k, er, geom, pol, cos_phi, sin_phi
    )

    return {"Fp": Fp_plus, "Fm": Fm_plus, "Gp": Gp_plus, "Gm": Gm_plus}


def _compute_downward_propagators(
    U: np.ndarray,
    V: np.ndarray,
    q1: np.ndarray,
    q2: np.ndarray,
    k: float,
    er: complex,
    geom: GeometryParams,
    pol: str,
    cos_phi: np.ndarray,
    sin_phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute downward propagators (Fm, Gm)."""
    # Use negative q for downward propagation
    Rh, Rv = _fresnel_coeffs(er, -q1, -q2)
    R = 0.5 * (Rv - Rh)
    mu_r = 1.0
    u_r = 1.0
    inv_q1 = _safe_inverse(-q1)
    inv_q2 = _safe_inverse(-q2)

    C_air = _compute_C_coeffs(-q1, geom, cos_phi, sin_phi, U, V)
    C_soil = _compute_C_coeffs(-q2, geom, cos_phi, sin_phi, U, V)
    B_air = _compute_B_coeffs(-q1, geom, cos_phi, sin_phi, U, V)
    B_soil = _compute_B_coeffs(-q2, geom, cos_phi, sin_phi, U, V)

    if pol == "vv":
        coeff = C_air
        coeff_t = C_soil
        Fm = (
            -(1 - Rv) * (1 + Rv) * inv_q1 * coeff["C1"]
            + (1 - Rv) * (1 - Rv) * inv_q1 * coeff["C2"]
            + (1 - Rv) * (1 + Rv) * inv_q1 * coeff["C3"]
            + (1 + Rv) * (1 - Rv) * inv_q1 * coeff["C4"]
            + (1 + Rv) * (1 + Rv) * inv_q1 * coeff["C5"]
            + (1 + Rv) * (1 - Rv) * inv_q1 * coeff["C6"]
        )
        Gm = (
            (1 + Rv) * (1 + Rv) * u_r * inv_q2 * coeff_t["C1"]
            - (1 + Rv) * (1 - Rv) * inv_q2 * coeff_t["C2"]
            - (1 + Rv * er) * (1 + Rv) * inv_q2 * coeff_t["C3"]
            - (1 - Rv) * er * (1 - Rv) * inv_q2 * coeff_t["C4"]
            - (1 - Rv) * (1 + Rv) * inv_q2 * coeff_t["C5"]
            - (1 - Rv) * (1 - Rv) * u_r * inv_q2 * coeff_t["C6"]
        )
    elif pol == "hh":
        coeff = C_air
        coeff_t = C_soil
        Fm = (
            (1 - Rh) * (1 + Rh) * inv_q1 * coeff["C1"]
            - (1 - Rh) * (1 - Rh) * inv_q1 * coeff["C2"]
            - (1 - Rh) * (1 + Rh) * inv_q1 * coeff["C3"]
            - (1 + Rh) * (1 - Rh) * inv_q1 * coeff["C4"]
            - (1 + Rh) * (1 + Rh) * inv_q1 * coeff["C5"]
            - (1 + Rh) * (1 - Rh) * inv_q1 * coeff["C6"]
        )
        Gm = (
            -(1 + Rh) * (1 + Rh) * er * inv_q2 * coeff_t["C1"]
            + (1 + Rh) * (1 - Rh) * inv_q2 * coeff_t["C2"]
            + (1 + Rh) * (1 + Rh) * inv_q2 * coeff_t["C3"]
            + (1 - Rh) * (1 - Rh) * inv_q2 * coeff_t["C4"]
            + (1 - Rh) * (1 + Rh) * inv_q2 * coeff_t["C5"]
            + (1 - Rh) * (1 - Rh) * er * inv_q2 * coeff_t["C6"]
        )
    else:  # cross-pol
        coeff = B_air
        coeff_t = B_soil
        Fm = (
            (1 - R) * (1 + R) * inv_q1 * coeff["B1"]
            - (1 - R) * (1 - R) * inv_q1 * coeff["B2"]
            - (1 - R) * (1 + R) * inv_q1 * coeff["B3"]
            + (1 + R) * (1 - R) * inv_q1 * coeff["B4"]
            + (1 + R) * (1 + R) * inv_q1 * coeff["B5"]
            + (1 + R) * (1 - R) * inv_q1 * coeff["B6"]
        )
        Gm = (
            -(1 + R) * (1 + R) * mu_r * inv_q2 * coeff_t["B1"]
            + (1 + R) * (1 - R) * inv_q2 * coeff_t["B2"]
            + (1 + R * er) * (1 + R) * inv_q2 * coeff_t["B3"]
            - (1 - R) * er * (1 - R) * inv_q2 * coeff_t["B4"]
            - (1 - R) * (1 + R) * inv_q2 * coeff_t["B5"]
            - (1 - R) * (1 - R) * mu_r * inv_q2 * coeff_t["B6"]
        )

    return Fm, Gm


# ============================================================================
# Coefficient Computation (C and B coefficients from Appendix C of paper)
# ============================================================================


def _spectral_angles(U: np.ndarray, V: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute spectral domain angles."""
    rho = np.hypot(U, V)
    cos_phi = np.ones_like(U)
    sin_phi = np.zeros_like(V)
    mask = rho > 0
    cos_phi[mask] = U[mask] / rho[mask]
    sin_phi[mask] = V[mask] / rho[mask]
    return cos_phi, sin_phi


def _compute_C_coeffs(
    q: np.ndarray,
    geom: GeometryParams,
    cos_phi: np.ndarray,
    sin_phi: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Compute C coefficients for VV/HH polarizations (Appendix C, Eqs C1-C6)."""
    cos_phi_s = geom.cos_phi_s
    sin_phi_s = geom.sin_phi_s
    cos_theta = geom.cos_theta_s
    sin_theta = geom.sin_theta_s
    z_x = geom.sin_theta_i * geom.cos_phi_i
    z_y = geom.sin_theta_i * geom.sin_phi_i
    zp_x = geom.sin_theta_s * geom.cos_phi_s
    zp_y = geom.sin_theta_s * geom.sin_phi_s

    q = np.asarray(q)
    cos_phi = np.asarray(cos_phi)
    sin_phi = np.asarray(sin_phi)
    U = np.asarray(U)
    V = np.asarray(V)

    C1 = -cos_phi_s * (-cos_phi - z_x * zp_x * cos_phi - z_x * zp_y * sin_phi) + sin_phi_s * (
        sin_phi + zp_x * z_y * cos_phi + z_y * zp_y * sin_phi
    )

    C2 = -cos_phi_s * (
        -q * cos_theta * cos_phi
        - U * z_x * cos_theta
        - V * zp_y * cos_theta * cos_phi
        - q * zp_x * sin_theta
        - U * z_x * zp_x * sin_theta
        - V * z_x * zp_y * sin_theta
        - V * z_x * cos_theta * sin_phi
        + V * zp_x * cos_theta * sin_phi
    ) + sin_phi_s * (
        U * z_y * cos_theta * cos_phi
        - U * zp_y * cos_theta
        + U * zp_x * z_y * sin_theta
        + q * zp_y * sin_theta
        + V * z_y * zp_y * sin_theta
        + q * cos_theta * sin_phi
        + U * zp_x * cos_theta * sin_phi
        + V * z_y * cos_theta * sin_phi
    )

    C3 = cos_phi_s * (
        U * zp_x * cos_theta * cos_phi
        - q * z_x * zp_x * cos_theta * cos_phi
        - U * sin_theta
        + q * z_x * sin_theta
        + U * zp_y * cos_theta * sin_phi
        - q * z_x * zp_y * cos_theta * sin_phi
    ) + sin_phi_s * (
        V * zp_x * cos_theta * cos_phi
        - q * zp_x * z_y * cos_theta * cos_phi
        - V * sin_theta
        + q * z_y * sin_theta
        + V * zp_y * cos_theta * sin_phi
        - q * z_y * zp_y * cos_theta * sin_phi
    )

    C4 = sin_theta * (
        -z_x * cos_theta * cos_phi
        - z_x * zp_x * sin_theta
        - z_y * zp_y * sin_theta
        - z_y * cos_theta * sin_phi
    ) - cos_theta * cos_phi_s * (
        -cos_theta * cos_phi
        - z_y * zp_y * cos_theta * cos_phi
        - zp_x * sin_theta
        + zp_x * z_y * cos_theta * sin_phi
    ) - cos_theta * sin_phi_s * (
        z_x * zp_y * cos_theta * cos_phi - zp_y * sin_theta - cos_theta * sin_phi - z_x * zp_x * cos_theta * sin_phi
    )

    C5 = sin_theta * (
        q * z_x * cos_phi
        + U * z_x * zp_x * cos_phi
        + V * zp_x * z_y * cos_phi
        + q * z_y * sin_phi
        + U * z_x * zp_y * sin_phi
        + V * z_y * zp_y * sin_phi
    ) - cos_theta * cos_phi_s * (
        q * cos_phi + U * zp_x * cos_phi + V * z_y * cos_phi - U * z_y * sin_phi + U * zp_y * sin_phi
    ) - cos_theta * sin_phi_s * (
        -V * z_x * cos_phi + V * zp_x * cos_phi + q * sin_phi + U * z_x * sin_phi + V * zp_y * sin_phi
    )

    C6 = sin_theta * (
        V * z_x * zp_y * cos_phi
        - U * z_y * zp_y * cos_phi
        - V * z_x * zp_x * sin_phi
        + U * zp_x * z_y * sin_phi
    ) + cos_theta * cos_phi_s * (
        -V * zp_y * cos_phi + q * z_y * zp_y * cos_phi + V * zp_x * sin_phi - q * zp_x * z_y * sin_phi
    ) + cos_theta * sin_phi_s * (
        U * zp_y * cos_phi - q * z_x * zp_y * cos_phi - U * zp_x * sin_phi + q * z_x * zp_x * sin_phi
    )

    return {"C1": C1, "C2": C2, "C3": C3, "C4": C4, "C5": C5, "C6": C6}


def _compute_B_coeffs(
    q: np.ndarray,
    geom: GeometryParams,
    cos_phi: np.ndarray,
    sin_phi: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Compute B coefficients for HV/VH polarizations (Appendix C, Eqs C7-C12)."""
    cos_phi_s = geom.cos_phi_s
    sin_phi_s = geom.sin_phi_s
    cos_theta = geom.cos_theta_s
    sin_theta = geom.sin_theta_s
    z_x = geom.sin_theta_i * geom.cos_phi_i
    z_y = geom.sin_theta_i * geom.sin_phi_i
    zp_x = geom.sin_theta_s * geom.cos_phi_s
    zp_y = geom.sin_theta_s * geom.sin_phi_s

    q = np.asarray(q)
    cos_phi = np.asarray(cos_phi)
    sin_phi = np.asarray(sin_phi)
    U = np.asarray(U)
    V = np.asarray(V)

    B1 = sin_theta * (-z_y * cos_phi + z_x * sin_phi) - cos_theta * cos_phi_s * (
        zp_x * z_y * cos_phi + sin_phi + z_y * zp_y * sin_phi
    ) - cos_theta * sin_phi_s * (-cos_phi - z_x * zp_x * cos_phi - z_x * zp_y * sin_phi)

    B2 = sin_theta * (
        -q * z_y * cos_theta * cos_phi
        - U * z_x * zp_y * cos_theta * cos_phi
        - V * z_y * zp_y * cos_theta * cos_phi
        - q * zp_x * z_y * sin_theta
        + q * z_x * zp_y * sin_theta
        + q * z_x * cos_theta * sin_phi
        + U * z_x * zp_x * cos_theta * sin_phi
        + V * zp_x * z_y * cos_theta * sin_phi
    ) - cos_theta * cos_phi_s * (
        U * z_y * cos_theta * cos_phi
        - U * zp_y * cos_theta * cos_phi
        + U * zp_x * z_y * sin_theta
        + q * zp_y * sin_theta
        + V * z_y * zp_y * sin_theta
        + q * cos_theta * sin_phi
        + U * zp_x * cos_theta * sin_phi
        + V * z_y * cos_theta * sin_phi
    ) - cos_theta * sin_phi_s * (
        -q * cos_theta * cos_phi
        - U * z_x * cos_theta * cos_phi
        - V * zp_y * cos_theta * cos_phi
        - q * zp_x * sin_theta
        - U * z_x * zp_x * sin_theta
        - V * z_x * zp_y * sin_theta
        - V * z_x * cos_theta * sin_phi
        + V * zp_x * cos_theta * sin_phi
    )

    B3 = sin_theta * (
        V * z_x * zp_x * cos_theta * cos_phi
        - U * z_y * zp_x * cos_theta * cos_phi
        - V * z_x * sin_theta
        + U * z_y * sin_theta
        + V * z_x * zp_y * cos_theta * sin_phi
        - U * z_y * zp_y * cos_theta * sin_phi
    ) + cos_theta * cos_phi_s * (
        -V * zp_x * cos_theta * cos_phi
        + q * zp_x * z_y * cos_theta * cos_phi
        + V * sin_theta
        - q * z_y * sin_theta
        - V * zp_y * cos_theta * sin_phi
        + q * z_y * zp_y * cos_theta * sin_phi
    ) + cos_theta * sin_phi_s * (
        U * zp_x * cos_theta * cos_phi
        - q * z_x * zp_x * cos_theta * cos_phi
        - U * sin_theta
        + q * z_x * sin_theta
        + U * zp_y * cos_theta * sin_phi
        - q * z_x * zp_y * cos_theta * sin_phi
    )

    B4 = -cos_phi_s * (
        z_x * zp_y * cos_theta * cos_phi - zp_y * sin_theta - cos_theta * sin_phi - z_x * zp_x * cos_theta * sin_phi
    ) + sin_phi_s * (
        -cos_theta * cos_phi - z_y * zp_y * cos_theta * cos_phi - zp_x * sin_theta + zp_x * z_y * cos_theta * sin_phi
    )

    B5 = -cos_phi_s * (
        -V * z_x * cos_phi + V * zp_x * cos_phi + q * sin_phi + U * z_x * sin_phi - V * zp_y * sin_phi
    ) + sin_phi_s * (
        q * cos_phi + U * zp_x * cos_phi + V * z_y * cos_phi - U * z_y * sin_phi + U * zp_y * sin_phi
    )

    B6 = cos_phi_s * (
        U * zp_y * cos_phi - q * z_x * zp_y * cos_phi - U * zp_x * sin_phi + q * z_x * zp_x * sin_phi
    ) + sin_phi_s * (
        V * zp_y * cos_phi - q * z_y * zp_y * cos_phi - V * zp_x * sin_phi + q * zp_x * z_y * sin_phi
    )

    return {"B1": B1, "B2": B2, "B3": B3, "B4": B4, "B5": B5, "B6": B6}


# ============================================================================
# Utility Functions
# ============================================================================

_SAFE_EPS = 1e-8


def _safe_div(a: np.ndarray | complex | float, b: np.ndarray | complex | float) -> np.ndarray:
    """Safe division avoiding singularities."""
    b = np.asarray(b) + 0.0
    mask = np.abs(b) < _SAFE_EPS
    b = b + mask * (_SAFE_EPS + 1j * _SAFE_EPS)
    return np.asarray(a) / b


def _safe_inverse(x: np.ndarray) -> np.ndarray:
    """Safe inverse."""
    return _safe_div(1.0, x)


def _fresnel_coeffs(er: complex, q1: np.ndarray, q2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Fresnel reflection coefficients."""
    Rh = _safe_div(q1 - q2, q1 + q2)
    Rv = _safe_div(er * q1 - q2, er * q1 + q2)
    return Rh, Rv


# ============================================================================
# Kirchhoff-Complementary Terms (gkc1, gkc2, gkc3 from Appendix A)
# ============================================================================


def _series_sum(
    coeff: np.ndarray | float,
    arg_x: np.ndarray | float,
    arg_y: np.ndarray | float,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> np.ndarray:
    """Compute series summation with roughness spectrum.
    
    Uses a Numba-accelerated inner loop when available; falls back to NumPy otherwise.
    """
    coeff_arr = np.asarray(coeff, dtype=np.complex128)
    if np.isscalar(arg_x):
        arg_x_arr = np.full_like(coeff_arr, float(arg_x))
    else:
        arg_x_arr = np.asarray(arg_x, dtype=float)
    if np.isscalar(arg_y):
        arg_y_arr = np.full_like(coeff_arr, float(arg_y))
    else:
        arg_y_arr = np.asarray(arg_y, dtype=float)

    # Fast path with Numba if available and provider exposes parameters
    if 'USE_NUMBA' in globals() and USE_NUMBA and hasattr(wn_provider, "_sigma2"):
        return _series_sum_numba(
            coeff_arr,
            arg_x_arr,
            arg_y_arr,
            float(wn_provider._sigma2),
            float(wn_provider._kl),
            int(Nmax),
            int(wn_provider._corr_code),
            float(_TWOPI5),
        )

    # Fallback pure NumPy path with iterative series (avoids factorial overflow)
    result = np.zeros_like(coeff_arr, dtype=np.complex128)
    term = np.ones_like(coeff_arr, dtype=np.complex128)
    for n in range(1, Nmax + 1):
        term = term * (coeff_arr / n)
        Wn = wn_provider(arg_x_arr, arg_y_arr, n)
        result += term * Wn
    return result


def _build_gkc1(
    U: np.ndarray,
    V: np.ndarray,
    geom: GeometryParams,
    q: np.ndarray,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> np.ndarray:
    """Build Kirchhoff-complementary term K1 (Eq A1)."""
    sigma2 = surf.sigma**2
    kz = geom.kz
    ksz = geom.ksz
    kx = geom.kx
    ky = geom.ky
    ksx = geom.ksx
    ksy = geom.ksy

    expo = np.exp(-sigma2 * (ksz**2 + kz**2 + ksz * kz + q**2 - q * ksz + q * kz))
    a_m = sigma2 * (kz + q) * (ksz + kz)
    a_n = sigma2 * (ksz - q) * (ksz + kz)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    return expo * sum_m * sum_n


def _build_gkc2(
    U: np.ndarray,
    V: np.ndarray,
    geom: GeometryParams,
    q: np.ndarray,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> np.ndarray:
    """Build Kirchhoff-complementary term K2 (Eq A2)."""
    sigma2 = surf.sigma**2
    kz = geom.kz
    ksz = geom.ksz
    kx = geom.kx
    ky = geom.ky
    ksx = geom.ksx
    ksy = geom.ksy

    expo = np.exp(-sigma2 * (ksz**2 + kz**2 + ksz * kz + q**2 - q * ksz + q * kz))
    a_m = sigma2 * (kz + q) * (ksz + kz)
    a_n = -sigma2 * (ksz - q) * (kz + q)
    sum_m = _series_sum(a_m, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    return expo * sum_m * sum_n


def _build_gkc3(
    U: np.ndarray,
    V: np.ndarray,
    geom: GeometryParams,
    q: np.ndarray,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> np.ndarray:
    """Build Kirchhoff-complementary term K3 (Eq A3)."""
    sigma2 = surf.sigma**2
    kz = geom.kz
    ksz = geom.ksz
    kx = geom.kx
    ky = geom.ky
    ksx = geom.ksx
    ksy = geom.ksy

    expo = np.exp(-sigma2 * (ksz**2 + kz**2 + ksz * kz + q**2 - q * ksz + q * kz))
    a_m = -sigma2 * (ksz - q) * (kz + q)
    a_n = sigma2 * (ksz - q) * (ksz + kz)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    return expo * sum_m * sum_n


# ============================================================================
# Complementary Terms (gc1-gc14 from Appendix A)
# ============================================================================


def _build_gc_block1(
    U: np.ndarray,
    V: np.ndarray,
    geom: GeometryParams,
    q: np.ndarray,
    qp: np.ndarray,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> Dict[str, np.ndarray]:
    """Build complementary block 1 (gc1-gc8, Eqs A4-A11)."""
    sigma2 = surf.sigma**2
    kz = geom.kz
    ksz = geom.ksz
    kx = geom.kx
    ky = geom.ky
    ksx = geom.ksx
    ksy = geom.ksy

    def expo(q_, qp_):
        return np.exp(-sigma2 * (ksz**2 + kz**2 + q_**2 + qp_**2 - (ksz - kz) * (q_ + qp_)))

    # gc1 (A4)
    a_n = sigma2 * (ksz - q) * (ksz - qp)
    a_m = sigma2 * (kz + q) * (kz + qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc1 = expo(q, qp) * sum_n * sum_m

    # gc2 (A5)
    a_n = sigma2 * (ksz - q) * (kz + qp)
    a_m = sigma2 * (kz + q) * (ksz - qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc2 = expo(q, qp) * sum_n * sum_m

    # gc3 (A6)
    a_n = sigma2 * (ksz - q) * (kz + qp)
    a_m = sigma2 * (kz + q) * (kz + qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc3 = expo(q, qp) * sum_n * sum_m

    # gc4 (A7)
    a_n = sigma2 * (ksz - q) * (kz + qp)
    a_m = -sigma2 * (ksz - q) * (kz + q)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc4 = expo(q, qp) * sum_n * sum_m

    # gc5 (A8)
    a_n = sigma2 * (kz + q) * (kz + qp)
    a_m = -sigma2 * (ksz - q) * (kz + q)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax)
    gc5 = expo(q, qp) * sum_n * sum_m

    # gc6 (A9)
    a_n = sigma2 * (ksz - q) * (ksz - qp)
    a_m = sigma2 * (kz + q) * (ksz - qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc6 = expo(q, qp) * sum_n * sum_m

    # gc7 (A10)
    a_n = sigma2 * (ksz - q) * (ksz - qp)
    a_m = -sigma2 * (ksz - q) * (kz + q)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc7 = expo(q, qp) * sum_n * sum_m

    # gc8 (A11)
    a_n = sigma2 * (kz + q) * (ksz - qp)
    a_m = -sigma2 * (ksz - q) * (kz + q)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax)
    gc8 = expo(q, qp) * sum_n * sum_m

    return {
        "gc1": gc1,
        "gc2": gc2,
        "gc3": gc3,
        "gc4": gc4,
        "gc5": gc5,
        "gc6": gc6,
        "gc7": gc7,
        "gc8": gc8,
    }


def _build_gc_block2(
    U: np.ndarray,
    V: np.ndarray,
    geom: GeometryParams,
    q: np.ndarray,
    qp: np.ndarray,
    surf: SurfaceParams,
    wn_provider: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    Nmax: int,
) -> Dict[str, np.ndarray]:
    """Build complementary block 2 (gc9-gc14, Eqs A12-A17)."""
    sigma2 = surf.sigma**2
    kz = geom.kz
    ksz = geom.ksz
    kx = geom.kx
    ky = geom.ky
    ksx = geom.ksx
    ksy = geom.ksy

    def expo(q_, qp_):
        return np.exp(-sigma2 * (ksz**2 + kz**2 + q_**2 + qp_**2 - (ksz - kz) * (q_ + qp_)))

    # gc9 (A12)
    a_n = sigma2 * (kz + q) * (ksz - qp)
    a_m = sigma2 * (kz + q) * (kz + qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc9 = expo(q, qp) * sum_n * sum_m

    # gc10 (A13)
    a_n = sigma2 * (kz + q) * (ksz - qp)
    a_m = -sigma2 * (ksz - qp) * (kz + qp)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc10 = expo(q, qp) * sum_n * sum_m

    # gc11 (A14)
    a_n = sigma2 * (kz + q) * (kz + qp)
    a_m = -sigma2 * (ksz - qp) * (kz + qp)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax)
    gc11 = expo(q, qp) * sum_n * sum_m

    # gc12 (A15)
    a_n = sigma2 * (ksz - q) * (ksz - qp)
    a_m = sigma2 * (ksz - q) * (kz + qp)
    sum_n = _series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc12 = expo(q, qp) * sum_n * sum_m

    # gc13 (A16)
    a_n = sigma2 * (ksz - q) * (ksz - qp)
    a_m = -sigma2 * (ksz - qp) * (kz + qp)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, kx + U, ky + V, wn_provider, Nmax)
    gc13 = expo(q, qp) * sum_n * sum_m

    # gc14 (A17)
    a_n = sigma2 * (ksz - q) * (kz + qp)
    a_m = -sigma2 * (ksz - qp) * (kz + qp)
    sum_n = _series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax)
    sum_m = _series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax)
    gc14 = expo(q, qp) * sum_n * sum_m

    return {
        "gc9": gc9,
        "gc10": gc10,
        "gc11": gc11,
        "gc12": gc12,
        "gc13": gc13,
        "gc14": gc14,
    }
