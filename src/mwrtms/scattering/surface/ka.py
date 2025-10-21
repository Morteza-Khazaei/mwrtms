"""Kirchhoff Approximation (KA) surface scattering model.

This module implements the Kirchhoff Approximation for rough surface bistatic
scattering with support for multiple autocorrelation functions (ACFs):
- Gaussian
- Exponential  
- x-Power (stretched-exponential)

The implementation follows the physical optics approach with a series expansion
in surface roughness. The only ACF-dependent component is the n-fold roughness
spectrum W^(n)(K).

References
----------
Based on the Kirchhoff/Physical Optics approximation for rough surface
scattering. See KA_MODEL_ACFs.md for complete derivation.
"""

from __future__ import annotations

import numpy as np

from ...factory import register_model
from ...core import PolarizationState, normalize_polarization
from .base import SurfaceScattering
from .iem.fresnel_utils import compute_fresnel_incident, compute_fresnel_specular
from .iem.geometry_utils import compute_spatial_frequency
from .iem.kirchhoff import VKA
from .iem.spectrum_aiem import compute_aiem_spectrum
from .iem.transition import compute_transition_function

__all__ = ["KAModel"]


@register_model("ka")
class KAModel(SurfaceScattering):
    """Kirchhoff Approximation with general autocorrelation functions.
    
    This model implements the complete bistatic KA scattering model with:
    - Gaussian, Exponential, or x-Power (stretched-exponential) ACFs
    - Vector Kirchhoff field coefficients (full polarimetric)
    - Series expansion in roughness moments with proper convergence
    - Exact Fresnel reflection coefficients
    
    The model is suitable for surfaces where:
    - kσ (normalized roughness) is moderate (typically < 3)
    - kL (normalized correlation length) is large (typically > 3)
    - Single scattering dominates (Kirchhoff approximation valid)
    
    Parameters
    ----------
    wave : ElectromagneticWave
        Electromagnetic wave properties
    geometry : ScatteringGeometry
        Scattering geometry (incident and scattered angles)
    surface : Surface
        Surface roughness properties (must provide rms_height and correlation_length)
    acf_type : str, optional
        Autocorrelation function type: 'gaussian', 'exponential', or 'xpower'
        (default: 'gaussian')
    alpha : float, optional
        Power exponent for x-power ACF (default: 1.5)
        Only used when acf_type='xpower'
    nmax : int, optional
        Optional cap on the series expansion. When ``None`` (default), the
        automatic AIEM convergence criterion selects the number of spectral
        terms.
    
    Notes
    -----
    The model computes the bistatic scattering coefficient as:
    
    σ⁰_qp = (k²/2) exp[-σ²(k_sz² + k_z²)] Σ(n=1 to nmax) (σ^(2n)/n!) |I_qp^(n)|² W^(n)(ΔK)
    
    where:
    - I_qp^(n) = q_z^n f_qp exp(-σ² |k_z| k_sz) is the KA moment kernel
    - W^(n)(K) is the n-fold roughness spectrum (ACF-dependent)
    - f_qp are the vector Kirchhoff field coefficients
    - q_z = k(cos θ_s + cos θ_i) is the vertical wavenumber sum
    
    ACF-specific spectra:
    - Gaussian: W^(n)(K) = (πL²/n) exp(-K²L²/4n)
    - Exponential: W^(n)(K) = 2πL²n / (n² + (KL)²)^(3/2)
    - x-Power: W^(n)(K) = L² n^(-2/α) Φ_α(KL n^(-1/α))
    
    References
    ----------
    Based on the Kirchhoff/Physical Optics approximation for rough surface
    scattering with Gaussian statistics. See KA_MODEL_ACFs.md for complete derivation.
    """

    MODEL_NAME = "KA"

    def __init__(
        self, 
        wave, 
        geometry, 
        surface=None, 
        *, 
        surface_roughness=None,
        acf_type: str = "exponential",
        alpha: float = 1.5,
        nmax: int | None = None
    ) -> None:
        """Initialize the KA model.
        
        Parameters
        ----------
        wave : ElectromagneticWave
            Electromagnetic wave properties
        geometry : ScatteringGeometry
            Scattering geometry
        surface : Surface, optional
            Surface roughness properties
        surface_roughness : Surface, optional
            Alternative parameter name for surface
        acf_type : str, optional
            Autocorrelation function: 'gaussian', 'exponential', or 'xpower'
        alpha : float, optional
            Power exponent for x-power ACF (default: 1.5)
        nmax : int, optional
            Optional cap on the AIEM spectral series order. If ``None`` (default),
            the automatic AIEM convergence criterion is used without additional
            truncation.
        """
        super().__init__(wave, geometry, surface, surface_roughness=surface_roughness)
        
        # Normalize ACF type
        self.acf_type = self._normalize_acf_type(acf_type)
        self.alpha = alpha
        self.nmax = nmax
        self._convergence_tolerance = 1e-16
        
        # Validate parameters
        if self._surface is None:
            raise ValueError("Surface roughness parameters required for KA model")
        
        if self.acf_type == "xpower" and alpha <= 0:
            raise ValueError("alpha must be positive for x-power ACF")
        
    @staticmethod
    def _normalize_acf_type(acf_type: str) -> str:
        """Normalize ACF type string."""
        acf_map = {
            'gaussian': 'gaussian',
            'gauss': 'gaussian',
            'exponential': 'exponential',
            'exp': 'exponential',
            'xpower': 'xpower',
            'x-power': 'xpower',
            'stretched': 'xpower',
            'stretched-exponential': 'xpower',
            'power': 'xpower',
            'alpha': 'xpower',
        }
        normalized = acf_map.get(acf_type.lower())
        if normalized is None:
            raise ValueError(f"Unknown ACF type: {acf_type}. "
                           f"Use 'gaussian', 'exponential', or 'xpower'")
        return normalized
        
    def compute(self, medium_above, medium_below, polarization) -> float:
        """Compute single-scattering Kirchhoff term aligned with AIEM."""
        pol_state = normalize_polarization(polarization)[0]
        
        # Surface statistics
        sigma = self._surface.rms_height()
        corr_length = self._surface.correlation_length()
        k = self._wave.wavenumber
        ks = k * sigma
        kl = k * corr_length
        
        # Geometry
        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_s_rad
        phi_i = 0.0
        phi_s = self._geometry.phi_s_rad
        cs = np.cos(theta_i)
        css = np.cos(theta_s)
        
        # Spectral order
        n_terms = self._determine_spectral_terms(ks, cs, css)
        if n_terms <= 0:
            return 0.0
        
        # Roughness spectrum (matches AIEM)
        spectra = self._compute_roughness_spectrum(
            kl, theta_i, theta_s, phi_s, phi_i, n_terms
        )
        
        # Fresnel reflection coefficients
        eps_r = medium_below.permittivity(self._wave.frequency_hz)
        Rvi, Rhi, _ = compute_fresnel_incident(eps_r, theta_i)
        Rvl, Rhl, _ = compute_fresnel_specular(eps_r, theta_i, theta_s, phi_s)
        
        # Transition-adjusted coefficients
        Tfv, Tfh = compute_transition_function(
            eps_r, theta_i, ks, cs, spectra, n_terms
        )
        Rvtran = Rvi + (Rvl - Rvi) * Tfv
        Rhtran = Rhi + (Rhl - Rhi) * Tfh
        
        # Kirchhoff field coefficients (vector KA)
        vka = VKA(theta_i, theta_s, phi_i, phi_s, Rvtran, Rhtran)
        fvv, fhh, fhv, fvh = vka.field_coefficients()
        
        sigma0 = self._compute_kirchhoff_term(
            fvv, fhh, fhv, fvh,
            ks, cs, css,
            spectra, n_terms,
            pol_state,
        )
        
        return float(max(sigma0, 0.0))

    def _compute_kirchhoff(self, R_h, R_v, polarization):
        """Return the Kirchhoff contribution (not used in this model)."""
        return 0.0

    def _compute_complementary(self, R_h, R_v, polarization):
        """Return the complementary contribution (not used in this model)."""
        return 0.0

    def _determine_spectral_terms(self, ks: float, cs: float, css: float) -> int:
        """Replicate AIEM spectral order selection for the Kirchhoff term."""
        if self.nmax is not None:
            if self.nmax <= 0:
                raise ValueError("nmax must be positive")
            max_terms = self.nmax
        else:
            max_terms = 1000
        
        cos_sum = cs + css
        factor = (ks ** 2) * (cos_sum ** 2)
        if factor <= 0.0:
            return 1
        
        tol = self._convergence_tolerance
        temp_prev = 0.0
        temp = factor
        iterm = 1
        
        while abs(temp - temp_prev) > tol and iterm < max_terms:
            temp_prev = temp
            iterm += 1
            temp = temp_prev * factor / iterm
        
        return max(1, min(iterm, max_terms))

    def _compute_roughness_spectrum(
        self,
        kl: float,
        theta_i: float,
        theta_s: float,
        phi_s: float,
        phi_i: float,
        n_terms: int,
    ) -> np.ndarray:
        """Compute AIEM roughness spectrum for the required orders."""
        if n_terms <= 0:
            return np.zeros(0, dtype=float)
        
        K = compute_spatial_frequency(kl, theta_i, theta_s, phi_s, phi_i)
        spectra = np.zeros(n_terms, dtype=float)
        corr_type = self._correlation_type_for_spectrum()
        
        for n in range(1, n_terms + 1):
            spectra[n - 1] = compute_aiem_spectrum(
                kl=kl,
                K=K,
                n=n,
                correlation_type=corr_type,
                power_exponent=self.alpha,
            )
        
        return spectra

    def _correlation_type_for_spectrum(self) -> str:
        """Map KA ACF identifiers to AIEM spectrum types."""
        if self.acf_type == "gaussian":
            return "gaussian"
        if self.acf_type == "exponential":
            return "exponential"
        if self.acf_type == "xpower":
            return "powerlaw"
        raise ValueError(f"Unsupported ACF type '{self.acf_type}' for spectrum computation.")

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
        polarization: PolarizationState,
    ) -> float:
        """Kirchhoff single-scattering term identical to AIEM implementation."""
        if n_terms <= 0 or len(spectra) < n_terms:
            return 0.0
        
        if polarization == PolarizationState.VV:
            field_coeff = fvv
        elif polarization == PolarizationState.HH:
            field_coeff = fhh
        elif polarization == PolarizationState.HV:
            field_coeff = fhv
        elif polarization == PolarizationState.VH:
            field_coeff = fvh
        else:
            return 0.0
        
        ks2 = ks ** 2
        cos_sum = cs + css
        sum_val = 0.0
        temp = 1.0
        
        for n in range(1, n_terms + 1):
            temp *= (ks2 * cos_sum ** 2) / n
            sum_val += temp * spectra[n - 1]
        
        exp_term = np.exp(-ks2 * cos_sum ** 2) * sum_val
        kterm = 0.5 * exp_term * np.abs(field_coeff) ** 2
        return float(np.real(kterm))
