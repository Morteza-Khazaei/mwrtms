"""Advanced Integral Equation Model (AIEM) for surface scattering.

This module implements the AIEM formulation with transition functions,
Kirchhoff terms, and complementary terms for all polarizations.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from .base import IEMBase, SurfaceRoughnessParameters
from .fresnel_utils import (
    compute_fresnel_incident,
    compute_fresnel_specular,
    compute_fresnel_nadir,
)
from .geometry_utils import compute_spatial_frequency
from .spectrum_aiem import compute_aiem_spectrum
from .transition import compute_transition_function
from .kirchhoff import compute_kirchhoff_coefficients, VKA
from .complementary import (
    compute_expal,
    compute_complementary_vv,
    compute_complementary_hh,
    compute_complementary_hv,
    compute_complementary_vh,
)
from .multiple_scattering import compute_multiple_scattering
from .guardrails import (
    validate_inputs,
    validate_wavenumbers,
    validate_fresnel_bounds,
    validate_field_coefficients,
    validate_sigma_real,
    validate_ms_balance,
    # Critical validations (Phase 1)
    validate_kirchhoff_polarization_independence,
    validate_cross_pol_kirchhoff_zero,
    validate_cross_pol_single_scattering,
    validate_cross_pol_magnitude,
    validate_cross_pol_ordering,
    validate_energy_conservation,
    ensure_reciprocity,
)
from ....core import PolarizationState
from ....medium import Medium
from ....factory import register_model

__all__ = ["AIEMModel"]


@register_model("aiem")
class AIEMModel(IEMBase):
    """Advanced Integral Equation Model (AIEM) for surface scattering.
    
    Implements the AIEM formulation with:
    - Transition function for reflection coefficients
    - Kirchhoff term for specular scattering
    - Complementary term for multiple scattering
    - Support for bistatic and backscatter geometries
    - All polarizations (VV, HH, HV, VH)
    
    Parameters
    ----------
    wave : ElectromagneticWave
        Electromagnetic wave properties
    geometry : ScatteringGeometry
        Scattering geometry (incident/scattered angles)
    surface : Surface
        Surface roughness model
    correlation_type : str, optional
        Surface correlation function ('gaussian', 'exponential', 'powerlaw')
        Default: 'exponential'
    correlation_length_m : float, optional
        Override correlation length from surface
    spectral_terms : int, optional
        Number of spectral terms (auto-computed if None)
    power_exponent : float, optional
        Exponent for power-law correlation (default 1.5)
    auto_terms : bool, optional
        Automatically determine spectral terms (default True)
    convergence_tolerance : float, optional
        Tolerance for series convergence (default 1e-16)
    include_multiple_scattering : bool, optional
        Include second-order multiple scattering (default False)
    ms_quadrature_points : int, optional
        Number of quadrature points for multiple scattering integration (default 129)
    ms_spectral_terms : int, optional
        Maximum spectral order for multiple scattering (default 8)
    enable_guardrails : bool, optional
        Activate physics and units guardrails for diagnostics (default True)
    
    References
    ----------
    Chen, K. S., Wu, T. D., Tsang, L., Li, Q., Shi, J., & Fung, A. K. (2003).
    Emission of rough surfaces calculated by the integral equation method with
    comparison to three-dimensional moment method simulations.
    IEEE Transactions on Geoscience and Remote Sensing, 41(1), 90-101.
    
    Wu, T. D., & Fung, A. K. (1992). A transition model for the reflection
    coefficient in surface scattering. IEEE TGRS, 30(4), 856-860.
    
    Examples
    --------
    >>> from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry
    >>> from mwrtms import build_surface_from_statistics, HomogeneousMedium
    >>> 
    >>> wave = ElectromagneticWave(5.4e9)  # 5.4 GHz
    >>> geometry = ScatteringGeometry(theta_i_deg=40.0)
    >>> surface = build_surface_from_statistics(
    ...     rms_height=0.01,  # 1 cm
    ...     correlation_length=0.05,  # 5 cm
    ...     correlation_type="exponential"
    ... )
    >>> 
    >>> model = AIEMModel(wave, geometry, surface)
    >>> air = HomogeneousMedium(1.0 + 0.0j)
    >>> soil = HomogeneousMedium(12.0 + 1.8j)
    >>> result = model.run(air, soil)
    >>> print(f"VV: {result.vv_db:.2f} dB")
    >>> print(f"HH: {result.hh_db:.2f} dB")
    """
    
    MODEL_NAME = "AIEM"
    
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
        convergence_tolerance: float = 1e-16,
        include_multiple_scattering: bool = False,
        ms_quadrature_points: int = 129,
        ms_spectral_terms: int = 8,
        ms_use_multiprocessing: bool = False,
        ms_workers: Optional[int] = None,
        enable_guardrails: bool = True,
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
        self._convergence_tolerance = convergence_tolerance
        self._include_multiple_scattering = include_multiple_scattering
        self._ms_quadrature_points = ms_quadrature_points
        self._ms_spectral_terms = ms_spectral_terms
        self._ms_use_multiprocessing = ms_use_multiprocessing
        self._ms_workers = ms_workers
        self._enable_guardrails = enable_guardrails
    
    def _compute_channel(
        self,
        medium_above: Medium,
        medium_below: Medium,
        polarization: PolarizationState,
        params: SurfaceRoughnessParameters,
    ) -> float:
        """Compute scattering coefficient for a single polarization channel."""
        
        # Extract parameters
        eps_r = medium_below.permittivity(self._wave.frequency_hz)
        k = self._wave.wavenumber
        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_s_rad
        phi_s = self._geometry.phi_s_rad
        phi_i = 0.0  # Standard convention
        
        ks = params.ks
        kl = params.kl
        sigma = params.sigma_m
        corr_length_m = params.correlation_length_m
        
        # Determine number of terms
        n_terms = self._determine_spectral_terms(params)
        
        # Compute roughness spectrum
        spectra = self._compute_roughness_spectrum(
            kl, theta_i, theta_s, phi_s, phi_i, n_terms
        )

        if self._enable_guardrails:
            validate_inputs(
                eps_r,
                self._wave.wavelength,
                sigma,
                corr_length_m,
                theta_i,
                theta_s,
            )
            validate_wavenumbers(k, theta_i, theta_s)
        
        # Compute Fresnel coefficients
        Rvi, Rhi, Rvhi = compute_fresnel_incident(eps_r, theta_i)
        # print(f"Rvi: {Rvi}, Rhi: {Rhi}, Rvhi: {Rvhi}")
        Rvl, Rhl, Rvhl = compute_fresnel_specular(eps_r, theta_i, theta_s, phi_s)

        if self._enable_guardrails:
            validate_fresnel_bounds(Rvi, Rhi, eps_r)
            validate_fresnel_bounds(Rvl, Rhl, eps_r)
        # print(f"Rvl: {Rvl}, Rhl: {Rhl}, Rvhl: {Rvhl}")
        # Rv0, Rh0 = compute_fresnel_nadir(eps_r)
        
        # Compute transition function
        cs = np.cos(theta_i)
        css = np.cos(theta_s)
        Tfv, Tfh = compute_transition_function(
            eps_r, theta_i, ks, cs, spectra, n_terms
        )
        
        # Transition-adjusted reflection coefficients
        Rvtran = Rvi + (Rvl - Rvi) * Tfv
        Rhtran = Rhi + (Rhl - Rhi) * Tfh
        Rvhtran = (Rvtran - Rhtran) / 2.0
        
        # Compute Kirchhoff field coefficients
        # fvv, fhh, fhv, fvh = compute_kirchhoff_coefficients(
        #     Rvtran, Rhtran, k, theta_i, theta_s, phi_s, phi_i
        # )

        vka_obj = VKA(theta_i, theta_s, phi_i, phi_s, Rvtran, Rhtran)
        fvv, fhh, fhv, fvh = vka_obj.field_coefficients()

        if self._enable_guardrails:
            validate_field_coefficients(
                {"vv": fvv, "hh": fhh, "hv": fhv, "vh": fvh}
            )
            # CRITICAL: Kirchhoff must be polarization-independent (Yang et al. 2017, p. 4741)
            validate_kirchhoff_polarization_independence(fhh, fvv)
        
        # Compute Kirchhoff term (single scattering)
        kterm = self._compute_kirchhoff_term(
            fvv, fhh, fhv, fvh, ks, cs, css, spectra, n_terms, polarization
        )
        
        # Phase 2: Validate Kirchhoff cross-pol vanishes in backscatter
        if self._enable_guardrails:
            # Check if backscatter geometry (theta_i ≈ theta_s, phi_s ≈ π)
            is_backscatter = (
                abs(theta_i - theta_s) < 1e-6 and 
                abs(phi_s - math.pi) < 1e-6
            )
            if is_backscatter and polarization in {PolarizationState.HV, PolarizationState.VH}:
                validate_cross_pol_kirchhoff_zero(kterm, threshold_db=-80.0)
        
        # Compute complementary term (single scattering)
        cterm = self._compute_complementary_term(
            Rvi, Rhi, Rvhi, eps_r, k, ks, sigma,
            theta_i, theta_s, phi_s, phi_i,
            fvv, fhh, fhv, fvh,
            spectra, n_terms, polarization
        )

        # Single scattering coefficient (Kirchhoff + complementary)
        sigma0_single = cterm #kterm + cterm # Using only complementary term because kirchhoff term is is being used in complementary term calculation
        if self._enable_guardrails:
            validate_sigma_real("single-scattering", sigma0_single)

        # Add multiple scattering if requested
        if self._include_multiple_scattering:
            kirchhoff_coeffs = {
                "vv": fvv,
                "hh": fhh,
                "hv": fhv,
                "vh": fvh,
            }
            ms_contrib = self._compute_multiple_scattering(
                eps_r,
                k,
                ks,
                kl,
                sigma,
                theta_i,
                theta_s,
                phi_i,
                phi_s,
                polarization=polarization,
                corr_length_m=corr_length_m,
                kirchhoff_coeffs=kirchhoff_coeffs,
            )
            if self._enable_guardrails:
                validate_sigma_real("multiple-scattering", ms_contrib)
                validate_ms_balance(
                    sigma0_single,
                    ms_contrib,
                    ks,
                    polarization.name.lower(),
                )
                # Cross-pol specific validation
                if polarization in {PolarizationState.HV, PolarizationState.VH}:
                    # Validate that cross-pol comes primarily from multiple scattering
                    validate_cross_pol_single_scattering(
                        sigma0_single,
                        ms_contrib,
                        threshold_ratio=0.01,
                    )
                    # Validate cross-pol magnitude is physically reasonable
                    validate_cross_pol_magnitude(
                        sigma0_single + ms_contrib,  # Total HV
                        ks,
                        max_reasonable_db=10.0,
                    )
            sigma0 = sigma0_single + ms_contrib
        else:
            sigma0 = sigma0_single
        if self._enable_guardrails:
            validate_sigma_real("total", sigma0)
        return float(np.real(sigma0))
    
    def _determine_spectral_terms(self, params: SurfaceRoughnessParameters) -> int:
        """Determine number of spectral terms for series convergence."""
        if not self._auto_terms:
            return self._spectral_terms or len(params.spectral_weights)
        
        # Use convergence criterion from MATLAB
        ks = params.ks
        cs = np.cos(self._geometry.theta_i_rad)
        css = np.cos(self._geometry.theta_s_rad)
        
        iterm = 1
        temp_old = 0.0
        temp = ks**2 * (cs + css)**2
        
        while abs(temp - temp_old) > self._convergence_tolerance:
            temp_old = temp
            iterm += 1
            temp = temp_old * (ks**2 * (cs + css)**2) / iterm
            if iterm > 1000:  # Safety limit
                break
        
        return iterm
    
    def _compute_roughness_spectrum(
        self,
        kl: float,
        theta_i: float,
        theta_s: float,
        phi_s: float,
        phi_i: float,
        n_terms: int
    ) -> np.ndarray:
        """Compute AIEM roughness spectrum for all orders."""
        
        # Compute spatial frequency component K
        K = compute_spatial_frequency(kl, theta_i, theta_s, phi_s, phi_i)
        
        spectra = np.zeros(n_terms)
        for n in range(1, n_terms + 1):
            spectra[n-1] = compute_aiem_spectrum(
                kl, K, n, self._correlation_type, self._power_exponent
            )
        
        return spectra
    
    def _compute_kirchhoff_term(
        self,
        fvv: complex,
        fhh: complex,
        fhv: complex,
        fvh: complex,
        ks: float,
        cs: float,
        css: float,
        spectra: np.ndarray,
        n_terms: int,
        polarization: PolarizationState
    ) -> float:
        """Compute Kirchhoff (specular) scattering term."""
        
        # Select field coefficient based on polarization
        if polarization == PolarizationState.VV:
            f = fvv
        elif polarization == PolarizationState.HH:
            f = fhh
        elif polarization == PolarizationState.HV:
            f = fhv
        elif polarization == PolarizationState.VH:
            f = fvh
        else:
            return 0.0
        
        # Compute series sum
        sum_val = 0.0
        temp = 1.0
        ks2 = ks**2
        
        for n in range(1, n_terms + 1):
            temp *= (ks2 * (cs + css)**2) / n
            sum_val += temp * spectra[n-1]
        
        # Kirchhoff term
        expk = np.exp(-ks2 * (css + cs)**2) * sum_val
        kterm = 0.5 * expk * np.abs(f)**2
        
        return float(kterm)
    
    def _compute_complementary_term(
        self,
        Rvi: complex,
        Rhi: complex,
        Rvhi: complex,
        eps_r: complex,
        k: float,
        ks: float,
        sigma: float,
        theta_i: float,
        theta_s: float,
        phi_s: float,
        phi_i: float,
        fvv: complex,
        fhh: complex,
        fhv: complex,
        fvh: complex,
        spectra: np.ndarray,
        n_terms: int,
        polarization: PolarizationState
    ) -> float:
        """Compute complementary (multiple scattering) term."""
        
        # Trigonometric quantities
        si = np.sin(theta_i)
        sis = np.sin(theta_s)
        cs = np.cos(theta_i)
        css = np.cos(theta_s)
        sfs = np.sin(phi_s)
        csfs = np.cos(phi_s)
        
        # Wave vector components
        # NOTE: These are magnitudes; sign/direction handled by direction parameter
        # in complementary field coefficient functions
        qq = cs
        qqt = np.sqrt(eps_r - si**2)
        qqs = css
        qqts = np.sqrt(eps_r - sis**2)
        
        # Compute scattering integrals based on polarization
        if polarization == PolarizationState.VV:
            I_pol = self._compute_scattering_integral_vv(
                Rvi, Rvhi, eps_r, k, ks, sigma,
                si, sis, cs, css, sfs, csfs,
                qq, qqt, qqs, qqts,
                fvv, spectra, n_terms
            )
        elif polarization == PolarizationState.HH:
            I_pol = self._compute_scattering_integral_hh(
                Rhi, Rvhi, eps_r, k, ks, sigma,
                si, sis, cs, css, sfs, csfs,
                qq, qqt, qqs, qqts,
                fhh, spectra, n_terms
            )
        elif polarization == PolarizationState.HV:
            I_pol = self._compute_scattering_integral_hv(
                Rvhi, eps_r, k, ks, sigma,
                si, sis, cs, css, sfs, csfs,
                qq, qqt, qqs, qqts,
                fhv, spectra, n_terms
            )
        elif polarization == PolarizationState.VH:
            I_pol = self._compute_scattering_integral_vh(
                Rvhi, eps_r, k, ks, sigma,
                si, sis, cs, css, sfs, csfs,
                qq, qqt, qqs, qqts,
                fvh, spectra, n_terms
            )
        else:
            return 0.0
        
        # Compute final complementary term
        sum_val = 0.0
        temp = 1.0
        ks2 = ks**2
        cs2 = cs**2
        css2 = css**2
        
        for n in range(1, n_terms + 1):
            temp *= (ks2 / n)
            sum_val += temp * np.abs(I_pol[n-1])**2 * spectra[n-1]
        
        cterm = 0.5 * np.exp(-ks2 * (cs2 + css2)) * sum_val
        
        return float(cterm)
    
    def _compute_scattering_integral_vv(
        self, Rv, Rvh, eps_r, k, ks, sigma,
        si, sis, cs, css, sfs, csfs,
        qq, qqt, qqs, qqts,
        fvv, spectra, n_terms
    ) -> np.ndarray:
        """Compute scattering integral for VV polarization."""
        
        Ivv = np.zeros(n_terms, dtype=np.complex128)
        ks2 = ks**2
        sigma2 = sigma**2
        
        # Wave vector components for complementary field
        qq1 = qq
        qq2 = qqs
        qq3 = qqt
        qq4 = qqts
        qq5 = qqt
        qq6 = qqts
        
        # Compute complementary field coefficients for all branches
        # Up/down for incident, up/down for scattered, air/substrate
        
        # Air-side, incident up
        Fvaupi = compute_complementary_vv(
            -si, 0.0, qq1, qq1, qq, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq1, ks, cs, css)
        
        # Air-side, incident down
        Fvadni = compute_complementary_vv(
            -si, 0.0, -qq1, -qq1, qq, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq1, ks, cs, css)
        
        # Air-side, scattered up
        Fvaups = compute_complementary_vv(
            -sis*csfs, -sis*sfs, qq2, qq2, qqs, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq2, ks, cs, css)
        
        # Air-side, scattered down
        Fvadns = compute_complementary_vv(
            -sis*csfs, -sis*sfs, -qq2, -qq2, qqs, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq2, ks, cs, css)
        
        # Substrate-side, incident up
        Fvbupi = compute_complementary_vv(
            -si, 0.0, qq3, qq5, qqt, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq5, ks, cs, css)
        
        # Substrate-side, incident down
        Fvbdni = compute_complementary_vv(
            -si, 0.0, -qq3, -qq5, qqt, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq5, ks, cs, css)
        
        # Substrate-side, scattered up
        Fvbups = compute_complementary_vv(
            -sis*csfs, -sis*sfs, qq4, qq6, qqts, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq6, ks, cs, css)
        
        # Substrate-side, scattered down
        Fvbdns = compute_complementary_vv(
            -sis*csfs, -sis*sfs, -qq4, -qq6, qqts, Rv, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq6, ks, cs, css)
        
        # Compute scattering integrals
        for n in range(1, n_terms + 1):
            fn = float(n)
            
            Ivv[n-1] = (
                (cs + css)**fn * fvv * np.exp(-ks2 * cs * css) +
                0.25 * (
                    Fvaupi * (css - qq1)**fn +
                    Fvadni * (css + qq1)**fn +
                    Fvaups * (cs + qq2)**fn +
                    Fvadns * (cs - qq2)**fn +
                    Fvbupi * (css - qq5)**fn +
                    Fvbdni * (css + qq5)**fn +
                    Fvbups * (cs + qq6)**fn +
                    Fvbdns * (cs - qq6)**fn
                )
            )
        
        return Ivv
    
    def _compute_scattering_integral_hh(
        self, Rh, Rvh, eps_r, k, ks, sigma,
        si, sis, cs, css, sfs, csfs,
        qq, qqt, qqs, qqts,
        fhh, spectra, n_terms
    ) -> np.ndarray:
        """Compute scattering integral for HH polarization."""
        
        Ihh = np.zeros(n_terms, dtype=np.complex128)
        ks2 = ks**2
        
        qq1 = qq
        qq2 = qqs
        qq5 = qqt
        qq6 = qqts
        
        # Compute complementary field coefficients (similar to VV but with Rh)
        Fhaupi = compute_complementary_hh(
            -si, 0.0, qq1, qq1, qq, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq1, ks, cs, css)
        
        Fhadni = compute_complementary_hh(
            -si, 0.0, -qq1, -qq1, qq, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq1, ks, cs, css)
        
        Fhaups = compute_complementary_hh(
            -sis*csfs, -sis*sfs, qq2, qq2, qqs, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq2, ks, cs, css)
        
        Fhadns = compute_complementary_hh(
            -sis*csfs, -sis*sfs, -qq2, -qq2, qqs, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq2, ks, cs, css)
        
        Fhbupi = compute_complementary_hh(
            -si, 0.0, qqt, qq5, qqt, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq5, ks, cs, css)
        
        Fhbdni = compute_complementary_hh(
            -si, 0.0, -qqt, -qq5, qqt, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq5, ks, cs, css)
        
        Fhbups = compute_complementary_hh(
            -sis*csfs, -sis*sfs, qqts, qq6, qqts, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq6, ks, cs, css)
        
        Fhbdns = compute_complementary_hh(
            -sis*csfs, -sis*sfs, -qqts, -qq6, qqts, Rh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq6, ks, cs, css)
        
        for n in range(1, n_terms + 1):
            Ihh[n-1] = (
                (cs + css)**n * fhh * np.exp(-ks2 * cs * css) +
                0.25 * (
                    Fhaupi * (css - qq1)**n +
                    Fhadni * (css + qq1)**n +
                    Fhaups * (cs + qq2)**n +
                    Fhadns * (cs - qq2)**n +
                    Fhbupi * (css - qq5)**n +
                    Fhbdni * (css + qq5)**n +
                    Fhbups * (cs + qq6)**n +
                    Fhbdns * (cs - qq6)**n
                )
            )
        
        return Ihh
    
    def _compute_scattering_integral_hv(
        self, Rvh, eps_r, k, ks, sigma,
        si, sis, cs, css, sfs, csfs,
        qq, qqt, qqs, qqts,
        fhv, spectra, n_terms
    ) -> np.ndarray:
        """Compute scattering integral for HV polarization."""
        
        Ihv = np.zeros(n_terms, dtype=np.complex128)
        ks2 = ks**2
        
        qq1 = qq
        qq2 = qqs
        qq5 = qqt
        qq6 = qqts
        
        # Compute complementary field coefficients for cross-pol
        Fhvaupi = compute_complementary_hv(
            -si, 0.0, qq1, qq1, qq, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq1, ks, cs, css)
        
        Fhvadni = compute_complementary_hv(
            -si, 0.0, -qq1, -qq1, qq, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq1, ks, cs, css)
        
        Fhvaups = compute_complementary_hv(
            -sis*csfs, -sis*sfs, qq2, qq2, qqs, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq2, ks, cs, css)
        
        Fhvadns = compute_complementary_hv(
            -sis*csfs, -sis*sfs, -qq2, -qq2, qqs, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq2, ks, cs, css)
        
        Fhvbupi = compute_complementary_hv(
            -si, 0.0, qqt, qq5, qqt, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq5, ks, cs, css)
        
        Fhvbdni = compute_complementary_hv(
            -si, 0.0, -qqt, -qq5, qqt, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq5, ks, cs, css)
        
        Fhvbups = compute_complementary_hv(
            -sis*csfs, -sis*sfs, qqts, qq6, qqts, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq6, ks, cs, css)
        
        Fhvbdns = compute_complementary_hv(
            -sis*csfs, -sis*sfs, -qqts, -qq6, qqts, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq6, ks, cs, css)
        
        for n in range(1, n_terms + 1):
            Ihv[n-1] = (
                (cs + css)**n * fhv * np.exp(-ks2 * cs * css) +
                0.25 * (
                    Fhvaupi * (css - qq1)**n +
                    Fhvadni * (css + qq1)**n +
                    Fhvaups * (cs + qq2)**n +
                    Fhvadns * (cs - qq2)**n +
                    Fhvbupi * (css - qq5)**n +
                    Fhvbdni * (css + qq5)**n +
                    Fhvbups * (cs + qq6)**n +
                    Fhvbdns * (cs - qq6)**n
                )
            )
        
        return Ihv
    
    def _compute_scattering_integral_vh(
        self, Rvh, eps_r, k, ks, sigma,
        si, sis, cs, css, sfs, csfs,
        qq, qqt, qqs, qqts,
        fvh, spectra, n_terms
    ) -> np.ndarray:
        """Compute scattering integral for VH polarization."""
        
        Ivh = np.zeros(n_terms, dtype=np.complex128)
        ks2 = ks**2
        
        qq1 = qq
        qq2 = qqs
        qq5 = qqt
        qq6 = qqts
        
        # Compute complementary field coefficients for cross-pol
        Fvhaupi = compute_complementary_vh(
            -si, 0.0, qq1, qq1, qq, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq1, ks, cs, css)
        
        Fvhadni = compute_complementary_vh(
            -si, 0.0, -qq1, -qq1, qq, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq1, ks, cs, css)
        
        Fvhaups = compute_complementary_vh(
            -sis*csfs, -sis*sfs, qq2, qq2, qqs, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(qq2, ks, cs, css)
        
        Fvhadns = compute_complementary_vh(
            -sis*csfs, -sis*sfs, -qq2, -qq2, qqs, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=False
        ) * compute_expal(-qq2, ks, cs, css)
        
        Fvhbupi = compute_complementary_vh(
            -si, 0.0, qqt, qq5, qqt, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq5, ks, cs, css)
        
        Fvhbdni = compute_complementary_vh(
            -si, 0.0, -qqt, -qq5, qqt, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq5, ks, cs, css)
        
        Fvhbups = compute_complementary_vh(
            -sis*csfs, -sis*sfs, qqts, qq6, qqts, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(qq6, ks, cs, css)
        
        Fvhbdns = compute_complementary_vh(
            -sis*csfs, -sis*sfs, -qqts, -qq6, qqts, Rvh, eps_r,
            si, sis, cs, css, sfs, csfs,
            is_substrate=True
        ) * compute_expal(-qq6, ks, cs, css)
        
        for n in range(1, n_terms + 1):
            Ivh[n-1] = (
                (cs + css)**n * fvh * np.exp(-ks2 * cs * css) +
                0.25 * (
                    Fvhaupi * (css - qq1)**n +
                    Fvhadni * (css + qq1)**n +
                    Fvhaups * (cs + qq2)**n +
                    Fvhadns * (cs - qq2)**n +
                    Fvhbupi * (css - qq5)**n +
                    Fvhbdni * (css + qq5)**n +
                    Fvhbups * (cs + qq6)**n +
                    Fvhbdns * (cs - qq6)**n
                )
            )
        
        return Ivh
    
    def _compute_multiple_scattering(
        self,
        eps_r: complex,
        k: float,
        ks: float,
        kl: float,
        sigma: float,
        theta_i: float,
        theta_s: float,
        phi_i: float,
        phi_s: float,
        polarization: PolarizationState,
        corr_length_m: float,
        kirchhoff_coeffs: dict[str, complex],
    ) -> float:
        """Compute second-order multiple scattering contribution.
        
        Parameters
        ----------
        eps_r : complex
            Relative permittivity
        k : float
            Wavenumber
        ks : float
            Normalized rms height
        kl : float
            Normalized correlation length
        sigma : float
            RMS height in meters
        corr_length_m : float
            Correlation length in meters
        theta_i : float
            Incident angle in radians
        theta_s : float
            Scattered angle in radians
        phi_i : float
            Incident azimuth in radians
        phi_s : float
            Scattered azimuth in radians
        polarization : PolarizationState
            Polarization state
        kirchhoff_coeffs : dict[str, complex]
            Kirchhoff field coefficients for each polarization
            
        Returns
        -------
        float
            Multiple scattering contribution (linear power)
        """
        # Map polarization to string
        pol_map = {
            PolarizationState.VV: "vv",
            PolarizationState.HH: "hh",
            PolarizationState.HV: "hv",
            PolarizationState.VH: "vh",
        }
        pol_str = pol_map.get(polarization, "vv")
        
        # Compute multiple scattering for this polarization
        # Note: Pass None for nmax to enable auto-determination based on roughness
        ms_results = compute_multiple_scattering(
            theta_i=theta_i,
            theta_s=theta_s,
            phi_i=phi_i,
            phi_s=phi_s,
            er=eps_r,
            ks=ks,
            kl=kl,
            k=k,
            sigma=sigma,
            surface_label=self._correlation_type,
            polarisations=[pol_str],
            n_points=self._ms_quadrature_points,
            nmax=self._ms_spectral_terms if self._ms_spectral_terms != 8 else 8,
        )
        
        return ms_results.get(pol_str, 0.0)
