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

import math
import numpy as np

from ...factory import register_model
from .base import SurfaceScattering
from .iem.fresnel_utils import compute_fresnel_incident
from .iem.geometry_utils import compute_q_vectors
from .iem.kirchhoff import compute_kirchhoff_coefficients
from .iem.spectrum_aiem import compute_aiem_spectrum

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
        Maximum order in series expansion (default: 8)
        Higher values improve accuracy but increase computation time
    
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
        nmax: int = 8
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
            Maximum order in series expansion (default: 8)
        """
        super().__init__(wave, geometry, surface, surface_roughness=surface_roughness)
        
        # Normalize ACF type
        self.acf_type = self._normalize_acf_type(acf_type)
        self.alpha = alpha
        self.nmax = nmax
        
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
        """Compute the backscatter coefficient for the given polarization.
        
        Parameters
        ----------
        medium_above : Medium
            Medium above the surface (typically air)
        medium_below : Medium
            Medium below the surface (typically soil/water)
        polarization : PolarizationState
            Desired polarization (VV, HH, HV, or VH)
            
        Returns
        -------
        float
            Backscatter coefficient (linear units)
        """
        # Get surface parameters
        sigma = self._surface.rms_height()  # meters
        L = self._surface.correlation_length()  # meters
        
        # Get permittivity
        eps_r = medium_below.permittivity(self._wave.frequency_hz)
        
        # Get wavelength
        lambda0 = self._wave.wavelength  # meters
        
        # For backscatter, incident and scattered angles are the same
        theta_i = self._geometry.theta_i_rad
        theta_s = self._geometry.theta_i_rad  # Backscatter
        phi_i = 0.0  # Incident azimuth (reference)
        phi_s = np.pi  # Backscatter azimuth (180 degrees)
        
        # Compute scattering coefficients for all polarizations
        result = self._sigma0_ka(
            lambda0, theta_i, phi_i, theta_s, phi_s, 
            sigma, L, eps_r, self.nmax
        )
        
        # Return the requested polarization
        pol_value = polarization.value.lower()
        if pol_value in ("vv", "v"):
            return float(max(result['VV'], 0.0))
        elif pol_value in ("hh", "h"):
            return float(max(result['HH'], 0.0))
        elif pol_value == "hv":
            return float(max(result['HV'], 0.0))
        elif pol_value == "vh":
            return float(max(result['VH'], 0.0))
        else:
            return 0.0

    def _compute_kirchhoff(self, R_h, R_v, polarization):
        """Return the Kirchhoff contribution (not used in this model)."""
        return 0.0

    def _compute_complementary(self, R_h, R_v, polarization):
        """Return the complementary contribution (not used in this model)."""
        return 0.0

    
    
    def _sigma0_ka(
        self,
        lambda0: float,
        theta_i: float,
        phi_i: float,
        theta_s: float,
        phi_s: float,
        sigma: float,
        L: float,
        eps_r: complex,
        nmax: int
    ) -> dict:
        """Compute KA scattering coefficients with general ACF.
        
        This implements the complete KA series expansion for bistatic scattering
        from a rough surface with the specified autocorrelation function.
        
        Uses IEM family components for:
        - Fresnel coefficients (compute_fresnel_incident)
        - Wave vector components (compute_q_vectors)
        - Kirchhoff field coefficients (compute_kirchhoff_coefficients)
        - Roughness spectrum (compute_aiem_spectrum)
        
        Parameters
        ----------
        lambda0 : float
            Wavelength (meters)
        theta_i : float
            Incident elevation angle (radians)
        phi_i : float
            Incident azimuth angle (radians)
        theta_s : float
            Scattered elevation angle (radians)
        phi_s : float
            Scattered azimuth angle (radians)
        sigma : float
            RMS height (meters)
        L : float
            Correlation length (meters)
        eps_r : complex
            Relative permittivity
        nmax : int
            Maximum series order
            
        Returns
        -------
        dict
            Dictionary with keys 'VV', 'HH', 'HV', 'VH' containing scattering
            coefficients in linear units
        """
        # Wavenumber
        k = 2.0 * np.pi / lambda0
        kl = k * L  # Normalized correlation length
        
        # Wave vector components using IEM geometry utilities
        kx, ky, ksx, ksy = compute_q_vectors(k, theta_i, theta_s, phi_s, phi_i)
        
        # Compute z-components manually (not provided by compute_q_vectors)
        kz = -k * np.cos(theta_i)  # Negative (downward)
        ksz = k * np.cos(theta_s)  # Positive (upward)
        
        # Horizontal wavenumber mismatch
        dKx = ksx - kx
        dKy = ksy - ky
        K = np.hypot(dKx, dKy)
        
        # Vertical wavenumber sum (crucial for KA)
        q_z = k * (np.cos(theta_s) + np.cos(theta_i))
        
        # Fresnel coefficients at incident angle using IEM utilities
        R_v, R_h, _ = compute_fresnel_incident(eps_r, theta_i)
        
        # Kirchhoff field coefficients using IEM utilities
        f_vv, f_hh, f_hv, f_vh = compute_kirchhoff_coefficients(
            R_v, R_h, k, theta_i, theta_s, phi_s, phi_i
        )
        
        # Outer exponential factor
        outer_exp = (k**2 / 2.0) * np.exp(-(sigma**2) * (ksz**2 + kz**2))
        
        # Inner exponential factor (in I_qp^(n))
        inner_exp = np.exp(-(sigma**2) * np.abs(kz) * ksz)
        
        # Series accumulation for each polarization
        def compute_series(f_qp: complex) -> float:
            """Compute series sum for given field coefficient."""
            accumulator = 0.0
            
            for n in range(1, nmax + 1):
                # KA moment kernel: I_qp^(n) = q_z^n * f_qp * exp(-σ² |k_z| k_sz)
                I_n = (q_z ** n) * f_qp * inner_exp
                
                # Roughness spectrum W^(n)(K) using IEM spectrum utilities
                # Map KA ACF types to AIEM correlation types
                if self.acf_type == "gaussian":
                    corr_type = "gaussian"
                elif self.acf_type == "exponential":
                    corr_type = "exponential"
                elif self.acf_type == "xpower":
                    corr_type = "powerlaw"
                else:
                    corr_type = "gaussian"
                
                W_n = compute_aiem_spectrum(
                    kl=kl,
                    K=K,
                    n=n,
                    correlation_type=corr_type,
                    power_exponent=self.alpha,
                    kx=kx,
                    ky=ky
                )
                
                # Series term: (σ^(2n) / n!) * |I_n|² * W^(n)
                factorial_n = math.factorial(n)
                term = (sigma ** (2 * n) / factorial_n) * (np.abs(I_n) ** 2) * W_n
                
                accumulator += term
            
            return float(np.real(outer_exp * accumulator))
        
        # Compute for all polarizations
        sigma0_vv = compute_series(f_vv)
        sigma0_hh = compute_series(f_hh)
        sigma0_hv = compute_series(f_hv)
        sigma0_vh = compute_series(f_vh)
        
        return {
            'VV': sigma0_vv,
            'HH': sigma0_hh,
            'HV': sigma0_hv,
            'VH': sigma0_vh
        }
