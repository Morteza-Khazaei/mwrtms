"""Kirchhoff field coefficients for AIEM."""

from __future__ import annotations

import numpy as np

__all__ = ["compute_kirchhoff_coefficients", "VKA"]


def compute_kirchhoff_coefficients(
    Rv: complex,
    Rh: complex,
    k: float,
    theta_i: float,
    theta_s: float,
    phi_s: float,
    phi_i: float = 0.0,
) -> tuple[complex, complex, complex, complex]:
    """
    Compute Kirchhoff field coefficients for all polarizations.
    
    Parameters
    ----------
    Rv : complex
        Vertical polarization reflection coefficient (transition-adjusted)
    Rh : complex
        Horizontal polarization reflection coefficient (transition-adjusted)
    k : float
        Wavenumber (2π/λ)
    theta_i : float
        Incident angle in radians
    theta_s : float
        Scattered angle in radians
    phi_s : float
        Scattered azimuth angle in radians
    phi_i : float, optional
        Incident azimuth angle in radians (default: 0.0)
    
    Returns
    -------
    fvv : complex
        VV polarization Kirchhoff field coefficient
    fhh : complex
        HH polarization Kirchhoff field coefficient
    fhv : complex
        HV polarization Kirchhoff field coefficient
    fvh : complex
        VH polarization Kirchhoff field coefficient
    
    Notes
    -----
    From AIEM.m, the Kirchhoff field coefficients are computed using
    the local surface slope components and reflection coefficients.
    
    The formulation accounts for:
    - Surface slope effects (zxx, zyy)
    - Polarization coupling
    - Geometric factors from incident and scattered directions
    """
    # Trigonometric quantities
    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    cs = np.cos(theta_i)
    css = np.cos(theta_s)
    cf = np.cos(phi_i)
    csfs = np.cos(phi_s)
    sfs = np.sin(phi_s)
    sf = np.sin(phi_i)
    
    cs2 = cs * cs
    si2 = si * si
    
    # Slope components
    denom = css + cs
    if np.abs(denom) < 1e-10:
        # Grazing angle - return zeros
        return 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j
    
    zxx = -(sis * csfs - si * cf) / denom
    zyy = -(sis * sfs - si * sf) / denom
    
    # Distance factor
    d2 = np.sqrt((zxx * cs - si * cf) ** 2 + (zyy * cs) ** 2 + (cs + si * zxx) ** 2)
    if np.abs(d2) < 1e-10:
        d2 = 1e-10
    
    # Field components for HH polarization
    hsnv = -(cs * csfs + si * (zxx * csfs + zyy * sfs))
    vsnh = css * csfs - zxx * sis
    hsnh = -sfs
    vsnv = zyy * cs * sis + css * (zyy * csfs * si - (cs + zxx * si) * sfs)
    
    # Tangential and normal components
    hsnt = (-(cs2 + si2) * sfs * (-si * cf + cs * zxx) + 
            csfs * (cs + si * zxx) * zyy + si * sfs * (zyy ** 2)) / d2
    
    hsnd = (-(cs + si * zxx) * (-csfs * si * cf + cs * csfs * zxx + cs * sfs * zyy)) / d2
    
    vsnt = ((cs2 + si2) * (-si * cf + cs * zxx) * (csfs * css - sis * zxx) + 
            css * sfs * (cs + si * zxx) * zyy - 
            (csfs * css * si * cf + cs * sis) * (zyy ** 2)) / d2
    
    vsnd = (-(cs + si * zxx) * (si * sis * zyy - 
            css * (si * sfs - cs * sfs * zxx + cs * csfs * zyy))) / d2
    
    # Kirchhoff field coefficients
    fhh = ((1.0 - Rh) * hsnv + (1.0 + Rh) * vsnh - 
           (hsnt + vsnd) * (Rh + Rv) * (zyy / d2))
    
    fvv = -((1.0 - Rv) * hsnv + (1.0 + Rv) * vsnh) + \
          (hsnt + vsnd) * (Rh + Rv) * (zyy / d2)
    
    fhv = (-(1.0 + Rv) * hsnh + (1.0 - Rv) * vsnv + 
           (hsnd - vsnt) * (Rh + Rv) * (zyy / d2))
    
    fvh = (-(1.0 + Rh) * hsnh + (1.0 - Rh) * vsnv + 
           (hsnd - vsnt) * (Rh + Rv) * (zyy / d2))
    
    return fvv, fhh, fhv, fvh


class VKA:
    """
    Vectorized Kirchhoff Approximation Model using NumPy vector operations.

    This implementation eliminates sign convention errors by using NumPy's
    built-in dot and cross product functions directly on vector arrays.
    """

    def __init__(self, theta_i, theta_s, phi_i, phi_s, Rv, Rh):
        """
        Initialize vectorized KA model.

        Parameters:
        -----------
        theta_i : float
            Incident elevation angle (radians)
        theta_s : float
            Scattered elevation angle (radians)
        phi_i : float
            Incident azimuth angle (radians)
        phi_s : float
            Scattered azimuth angle (radians)
        Rv : complex
            Vertical polarization Fresnel reflection coefficient
        Rh : complex
            Horizontal polarization Fresnel reflection coefficient
        """
        self.theta_i = theta_i
        self.theta_s = theta_s
        self.phi_i = phi_i
        self.phi_s = phi_s
        self.Rv = Rv
        self.Rh = Rh

        # Precompute all vectors using exact definitions
        self._setup_vectors()

        # Compute surface slopes and normalization factors
        self.zx, self.zy = self._compute_surface_slopes()
        self.D1 = np.sqrt(1 + self.zx**2 + self.zy**2)
        self.D2 = np.sqrt(self.zy**2 + (np.sin(theta_i) - self.zx * np.cos(theta_i))**2)

        # Update surface-dependent vectors
        self._setup_surface_vectors()

    def _setup_vectors(self):
        """Define all polarization and wave vectors using exact definitions."""

        # Wave vectors (following the paper's convention)
        # Incident: downward propagation (negative z-component)
        self.k_i = np.array([
            np.sin(self.theta_i) * np.cos(self.phi_i),
            np.sin(self.theta_i) * np.sin(self.phi_i),
            -np.cos(self.theta_i)  # Negative for downward propagation
        ])

        # Scattered: upward propagation (positive z-component)
        self.k_s = np.array([
            np.sin(self.theta_s) * np.cos(self.phi_s),
            np.sin(self.theta_s) * np.sin(self.phi_s),
            np.cos(self.theta_s)   # Positive for upward propagation
        ])

        # Incident polarization vectors (from paper definitions)
        self.h_i = np.array([
            -np.sin(self.phi_i),
            np.cos(self.phi_i),
            0
        ])

        # Use the paper's definition for v_i (with the negative sign)
        self.v_i = np.array([
            np.cos(self.theta_i) * np.cos(self.phi_i),
            np.cos(self.theta_i) * np.sin(self.phi_i),
            np.sin(self.theta_i)
        ])

        # Scattered polarization vectors (from paper definitions)
        self.h_s = -np.array([
            -np.sin(self.phi_s),
            np.cos(self.phi_s),
            0
        ])

        self.v_s = -np.array([
            np.cos(self.theta_s) * np.cos(self.phi_s),
            np.cos(self.theta_s) * np.sin(self.phi_s),
            -np.sin(self.theta_s)
        ])

    def _compute_surface_slopes(self):
        """Compute surface slopes using stationary phase approximation."""
        # Wave vector components
        kx_i = self.k_i[0]
        ky_i = self.k_i[1]
        kz_i = self.k_i[2]

        kx_s = self.k_s[0]
        ky_s = self.k_s[1]
        kz_s = self.k_s[2]

        # Stationary phase relations (Eq. 18)
        denominator = kz_s - kz_i
        if np.isclose(denominator, 0.0):
            raise ValueError("Invalid scattering geometry: kz_s - kz_i ≈ 0")

        # Stationary-phase slopes (Eq. 20): zx = -(k_sx - k_ix)/(k_sz - k_iz)
        zx = -(kx_s - kx_i) / denominator
        zy = -(ky_s - ky_i) / denominator

        return zx, zy

    def _setup_surface_vectors(self):
        """Define surface normal and local surface vectors."""

        # Surface normal vector
        self.n = np.array([-self.zx, -self.zy, 1]) / self.D1

        # Compute t = (k_i × n) / |k_i × n|
        k_cross_n = np.cross(self.k_i, self.n)
        k_cross_n_magnitude = np.linalg.norm(k_cross_n)
        self.t = k_cross_n / k_cross_n_magnitude

        # Compute d = k_i × t
        self.d = np.cross(self.k_i, self.t)

    def compute_all_dot_products(self):
        """
        Compute all required dot and cross products using NumPy operations.
        This eliminates any manual calculation errors and sign convention issues.
        """
        results = {}

        # Cross products with surface normal
        results['vs_dot_n_cross_vi'] = np.dot(self.v_s, np.cross(self.n, self.v_i))
        results['vs_dot_n_cross_hi'] = np.dot(self.v_s, np.cross(self.n, self.h_i))
        results['hs_dot_n_cross_vi'] = np.dot(self.h_s, np.cross(self.n, self.v_i))
        results['hs_dot_n_cross_hi'] = np.dot(self.h_s, np.cross(self.n, self.h_i))

        # Dot products with local surface vectors
        results['vi_dot_t'] = np.dot(self.v_i, self.t)
        results['hi_dot_d'] = np.dot(self.h_i, self.d)
        results['vs_dot_t'] = np.dot(self.v_s, self.t)
        results['vs_dot_d'] = np.dot(self.v_s, self.d)
        results['hs_dot_t'] = np.dot(self.h_s, self.t)
        results['hs_dot_d'] = np.dot(self.h_s, self.d)

        # Other required dot products
        results['n_dot_ki'] = np.dot(self.n, self.k_i)
        results['n_dot_d'] = np.dot(self.n, self.d)
        results['hs_dot_ki'] = np.dot(self.h_s, self.k_i)
        results['vs_dot_ki'] = np.dot(self.v_s, self.k_i)

        return results

    def field_coefficients(self):
        """
        Compute Kirchhoff field coefficients using vectorized operations.

        Returns:
        --------
        tuple : (fhh, fvv, fvh, fhv)
            Complex field coefficients
        """
        # Get all dot products
        dots = self.compute_all_dot_products()

        # Extract reflection coefficients
        Rv, Rh = self.Rv, self.Rh
        D1 = self.D1

        # co-pol common term
        co_pol_bracket = (dots['hs_dot_d'] * dots['n_dot_ki']
                          - dots['n_dot_d'] * dots['hs_dot_ki']
                          - dots['vs_dot_t'] * dots['n_dot_ki'])

        # x-pol common term
        xpol_bracket = (dots['hs_dot_t'] * dots['n_dot_ki']
                        - dots['n_dot_d'] * dots['vs_dot_ki']
                        + dots['vs_dot_d'] * dots['n_dot_ki'])

        # f_vv: Vertical incident → Vertical scattered
        fvv = -((1 - Rv) * dots['hs_dot_n_cross_vi'] + (1 + Rv) * dots['vs_dot_n_cross_hi']) * D1 \
            - (Rh + Rv) * dots['vi_dot_t'] * co_pol_bracket * D1

        # f_vh: Vertical incident → Horizontal scattered
        fvh = ((1 - Rh) * dots['vs_dot_n_cross_vi'] - (1 + Rh) * dots['hs_dot_n_cross_hi']) * D1 \
            - (Rh + Rv) * dots['hi_dot_d'] * xpol_bracket * D1

        # f_hv: Horizontal incident → Vertical scattered
        fhv = -((1 - Rv) * dots['vs_dot_n_cross_vi'] - (1 + Rv) * dots['hs_dot_n_cross_hi']) * D1 \
            -(Rh + Rv) * dots['vi_dot_t'] * xpol_bracket * D1

        # f_hh: Horizontal incident → Horizontal scattered
        fhh = ((1 + Rh) * dots['vs_dot_n_cross_hi'] + (1 - Rh) * dots['hs_dot_n_cross_vi']) * D1 \
            - (Rh + Rv) * dots['hi_dot_d'] * co_pol_bracket * D1

        return fhh, fvv, fvh, fhv

    def verify_orthogonality(self):
        """Verify electromagnetic orthogonality constraints using NumPy."""
        return {
            'hi_dot_ki': np.dot(self.h_i, self.k_i),
            'vi_dot_ki': np.dot(self.v_i, self.k_i),
            'hi_dot_vi': np.dot(self.h_i, self.v_i),
            'hs_dot_ks': np.dot(self.h_s, self.k_s),
            'vs_dot_ks': np.dot(self.v_s, self.k_s),
            'hs_dot_vs': np.dot(self.h_s, self.v_s),
            'norm_hi': np.linalg.norm(self.h_i),
            'norm_vi': np.linalg.norm(self.v_i),
            'norm_hs': np.linalg.norm(self.h_s),
            'norm_vs': np.linalg.norm(self.v_s)
        }

    def compare_with_analytical(self):
        """
        Compare vectorized results with our analytical derivations.
        This helps verify which analytical expressions were correct.
        """
        dots = self.compute_all_dot_products()

        print("=== VECTORIZED vs ANALYTICAL COMPARISON ===\n")

        # Compare key expressions we derived
        print("Cross products with surface normal:")
        print(f"v_s · (n × v_i) = {dots['vs_dot_n_cross_vi']:.6f}")
        print(f"v_s · (n × h_i) = {dots['vs_dot_n_cross_hi']:.6f}")
        print(f"h_s · (n × v_i) = {dots['hs_dot_n_cross_vi']:.6f}")
        print(f"h_s · (n × h_i) = {dots['hs_dot_n_cross_hi']:.6f}")

        print("\nDot products with local vectors:")
        print(f"v_i · t = {dots['vi_dot_t']:.6f}")
        print(f"h_i · d = {dots['hi_dot_d']:.6f}")
        print(f"v_s · t = {dots['vs_dot_t']:.6f}")
        print(f"v_s · d = {dots['vs_dot_d']:.6f}")
        print(f"h_s · t = {dots['hs_dot_t']:.6f}")
        print(f"h_s · d = {dots['hs_dot_d']:.6f}")

        print("\nOther dot products:")
        print(f"n · k_i = {dots['n_dot_ki']:.6f}")
        print(f"n · d = {dots['n_dot_d']:.6f}")
        print(f"h_s · k_i = {dots['hs_dot_ki']:.6f}")
        print(f"v_s · k_i = {dots['vs_dot_ki']:.6f}")

        return dots

    def get_all_vectors(self):
        """Return all computed vectors for inspection."""
        return {
            'k_i': self.k_i,
            'k_s': self.k_s,
            'h_i': self.h_i,
            'v_i': self.v_i,
            'h_s': self.h_s,
            'v_s': self.v_s,
            'n': self.n,
            't': self.t,
            'd': self.d,
            'surface_slopes': (self.zx, self.zy),
            'normalization': (self.D1, self.D2)
        }
