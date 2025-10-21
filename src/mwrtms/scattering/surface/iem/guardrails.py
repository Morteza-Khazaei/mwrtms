"""Enhanced physics and mathematical guardrails for AIEM multiple scattering.

This module provides comprehensive validation for AIEM implementations based on:
1. First-principles electromagnetic theory
2. Gaussian surface statistics
3. Perturbation theory convergence requirements
4. Energy conservation and reciprocity
5. Spectral coupling structure from Yang et al. (2017)

Critical additions based on equation re-derivation:
- Negative integrand values from destructive interference are PHYSICAL
- Propagator evaluation at shifted spectral points must be validated
- Phase factor structure must match Gaussian moment generating functions
- Series convergence requires careful term-by-term analysis
- Integration domain must capture interference side lobes (not just main peak)
"""

from __future__ import annotations

import math
import warnings
from typing import Mapping, Optional, Tuple, Dict

import numpy as np

from .kirchhoff import VKA
from .complementary import _compute_slope_and_field_components

__all__ = [
    "GuardrailViolation",
    # Basic validation
    "validate_inputs",
    "validate_wavenumbers",
    "validate_fresnel_bounds",
    "validate_field_coefficients",
    "ensure_array_finite",
    "validate_sigma_real",
    # Multiple scattering specific
    "validate_ms_balance",
    "validate_cross_pol_single_scattering",
    "validate_cross_pol_magnitude",
    "validate_cross_pol_ordering",
    "validate_kirchhoff_polarization_independence",
    "validate_cross_pol_kirchhoff_zero",
    "validate_component_hierarchy",
    "validate_hh_vv_ms_asymmetry",
    # NEW: Enhanced phase and propagator validation
    "validate_phase_factor_structure",
    "validate_propagator_shifts",
    "validate_spectral_coupling_structure",
    "validate_propagator_conjugation",
    # NEW: Enhanced integrand validation
    "validate_integrand_structure",
    "classify_negative_integrand",
    "validate_integration_domain",
    # NEW: Enhanced series convergence
    "validate_double_summation_structure",
    "validate_series_term_ratios",
    "detect_series_divergence",
    # Reciprocity and symmetry
    "ensure_reciprocity",
    "validate_nadir_symmetry",
    "validate_brewster_behavior",
    # Physical bounds
    "validate_energy_conservation",
    "validate_cross_pol_ordering",
    "validate_angular_trends",
    # Integration and convergence
    "assert_spectrum_consistency",
    "enforce_nonnegative_integrand",
    "enforce_min_integrand_energy",
    "validate_radiation_condition",
    "validate_series_convergence",
    "validate_double_summation_convergence",
    # Field coefficients
    "validate_field_coefficient_structure",
    "validate_complementary_coefficients",
    "validate_pec_limit",
    # Spectral coupling
    "validate_spectral_coupling",
    "validate_phase_factors",
    # Limiting cases
    "validate_spm_regime",
    "validate_smooth_limit",
    # Comprehensive validation suites
    "validate_aiem_single_scattering",
    "validate_aiem_multiple_scattering",
    "validate_aiem_multiple_scattering_comprehensive",
    "validate_aiem_total_scattering",
    # Diagnostics
    "diagnose_vka_complementary_slope_alignment",
]


class GuardrailViolation(RuntimeError):
    """Raised when a physics/mathematics guardrail is violated."""


def _guard(condition: bool, message: str, *, severity: str = "error") -> None:
    """Check condition and raise/warn if violated."""
    if condition:
        return
    if severity == "warning":
        warnings.warn(message, UserWarning, stacklevel=3)
    else:
        raise GuardrailViolation(message)


# =============================================================================
# SECTION 1: ENHANCED INTEGRAND VALIDATION
# =============================================================================


def validate_integrand_structure(
    integrand: np.ndarray,
    label: str,
    expected_peak_location: Optional[Tuple[float, float]] = None,
    ks: float = 0.5,
) -> None:
    """Validate that integrand has expected spectral structure.
    
    Based on re-derivation: Integrands should:
    1. Peak near specular direction (for Kirchhoff-dominated scattering)
    2. Have side lobes from interference (for multiple scattering)
    3. Decay at large |u|, |v| due to radiation condition
    
    Parameters
    ----------
    integrand : np.ndarray
        Complex or real integrand array (Nu, Nv)
    label : str
        Identifier for error messages
    expected_peak_location : Optional[Tuple[float, float]]
        Expected (u, v) location of peak (e.g., specular direction)
    ks : float
        Normalized roughness for context
    """
    integrand_real = np.real(integrand)
    integrand_abs = np.abs(integrand)
    
    # Check 1: Should have some non-zero energy
    max_abs = float(np.max(integrand_abs))
    _guard(
        max_abs > 1e-30,
        f"{label}: Integrand is essentially zero (max = {max_abs:.3e}). "
        "This indicates missing physics or numerical underflow.",
    )


# =============================================================================
# DIAGNOSTICS
# =============================================================================


def diagnose_vka_complementary_slope_alignment(
    theta_i: float,
    theta_s: float,
    phi_i: float,
    phi_s: float,
    *,
    tolerance: float = 1e-6,
    emit_warning: bool = True,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Compare VKA stationary-phase slopes against complementary-field slopes at (u, v) = (0, 0).

    This diagnostic evaluates the invariant that, for zero spectral offset, the complementary
    slope components should agree with the stationary-phase (VKA) slopes when the vertical
    wavenumber denominators are well-defined. Branches with nearly singular denominators are
    reported but excluded from the tolerance check.

    Parameters
    ----------
    theta_i, theta_s : float
        Incident and scattered elevation angles (radians).
    phi_i, phi_s : float
        Incident and scattered azimuth angles (radians).
    tolerance : float, optional
        Acceptable absolute tolerance for slope agreement when denominators are regular.
    emit_warning : bool, optional
        If True, emit a warning when no complementary branch agrees with the VKA slopes.

    Returns
    -------
    Dict[str, Dict[str, Dict[str, float]]]
        Summary containing the VKA slopes and per-branch complementary slopes with the
        resulting differences.
    """
    vka_model = VKA(theta_i, theta_s, phi_i, phi_s, Rv=0.0 + 0.0j, Rh=0.0 + 0.0j)
    zx_vka = float(np.real(vka_model.zx))
    zy_vka = float(np.real(vka_model.zy))

    si = np.sin(theta_i)
    sis = np.sin(theta_s)
    cs = np.cos(theta_i)
    css = np.cos(theta_s)
    sfs = np.sin(phi_s)
    csfs = np.cos(phi_s)

    branches = {
        "air_incident_up": (cs, cs),
        "air_incident_down": (-cs, -cs),
        "air_scattered_up": (css, css),
        "air_scattered_down": (-css, -css),
    }

    branch_results: Dict[str, Dict[str, float]] = {}
    aligned = False

    for label, (q, qslp) in branches.items():
        denom_scatter = css - qslp
        singular = abs(denom_scatter) < 1e-12

        if singular:
            branch_results[label] = {
                "singular": 1.0,
                "zx": float("nan"),
                "zy": float("nan"),
                "delta_norm": float("nan"),
            }
            continue

        zx_comp, zy_comp, *_ = _compute_slope_and_field_components(
            0.0, 0.0, q, qslp, si, sis, cs, css, sfs, csfs
        )
        zx_comp_real = float(np.real(zx_comp))
        zy_comp_real = float(np.real(zy_comp))

        delta_x = zx_comp_real - zx_vka
        delta_y = zy_comp_real - zy_vka
        delta_norm = float(math.hypot(delta_x, delta_y))

        if delta_norm <= tolerance:
            aligned = True

        branch_results[label] = {
            "singular": 0.0,
            "zx": zx_comp_real,
            "zy": zy_comp_real,
            "delta_norm": delta_norm,
        }

    if emit_warning and not aligned:
        warnings.warn(
            "No complementary slope branch matched the VKA stationary-phase slopes "
            f"within tolerance {tolerance:.2e}. Inspect branch diagnostics for details.",
            UserWarning,
            stacklevel=2,
        )

    return {
        "vka": {"zx": zx_vka, "zy": zy_vka},
        "branches": branch_results,
    }
    
    # Check 2: Should not be constant (would indicate missing spatial variation)
    std_abs = float(np.std(integrand_abs))
    _guard(
        std_abs / max_abs > 1e-6,
        f"{label}: Integrand has no spatial variation (std/max = {std_abs/max_abs:.3e}). "
        "This suggests missing spectral coupling or incorrect spectrum computation.",
        severity="warning",
    )
    
    # Check 3: For rough surfaces, should show interference structure
    if ks > 0.5:
        # Count zero crossings in real part (proxy for interference fringes)
        real_signs = np.sign(integrand_real)
        sign_changes = np.sum(np.abs(np.diff(real_signs, axis=0)) > 0)
        sign_changes += np.sum(np.abs(np.diff(real_signs, axis=1)) > 0)
        
        _guard(
            sign_changes > integrand.shape[0] // 4,
            f"{label}: Too few sign changes ({sign_changes}) for kσ = {ks:.2f}. "
            "Rough surfaces should show interference fringes in integrand. "
            "Check if negative series coefficients are being preserved.",
            severity="warning",
        )
    
    # Check 4: Should decay at edges (radiation condition)
    edge_max = max(
        float(np.max(integrand_abs[0, :])),
        float(np.max(integrand_abs[-1, :])),
        float(np.max(integrand_abs[:, 0])),
        float(np.max(integrand_abs[:, -1])),
    )
    center_max = float(np.max(integrand_abs[
        integrand.shape[0]//4:-integrand.shape[0]//4,
        integrand.shape[1]//4:-integrand.shape[1]//4
    ]))
    
    if center_max > 0:
        edge_ratio = edge_max / center_max
        _guard(
            edge_ratio < 0.5,
            f"{label}: Integrand does not decay at edges (edge/center = {edge_ratio:.3f}). "
            "Integration domain may be too small or radiation condition not applied properly.",
            severity="warning",
        )


def classify_negative_integrand(
    integrand: np.ndarray,
    label: str,
    reference_scale: float,
    ks: float,
) -> Dict[str, float]:
    """Classify negative integrand values as physical or numerical error.
    
    CRITICAL INSIGHT from re-derivation:
    Negative integrand values arise from TWO sources:
    
    1. PHYSICAL: Destructive interference from negative series coefficients
       - Equations A7, A8, A10, A11, A13, A14 have -σ²(...) terms
       - These create interference: 1 - σ²A + (σ²A)²/2! - ...
       - Expected magnitude: O(10^-6 to 10^-4) relative to peak
       - Sign pattern: Should be spatially coherent (not random noise)
    
    2. NUMERICAL ERROR: Accumulated floating-point errors
       - Expected magnitude: O(10^-12 to 10^-10) relative to peak  
       - Sign pattern: Random spatial distribution
    
    Parameters
    ----------
    integrand : np.ndarray
        Real integrand array
    label : str
        Identifier
    reference_scale : float
        Peak magnitude for normalization
    ks : float
        Normalized roughness
        
    Returns
    -------
    Dict[str, float]
        Classification metrics
    """
    integrand_real = np.real(integrand)
    negative_mask = integrand_real < 0
    
    if not np.any(negative_mask):
        return {
            "has_negatives": False,
            "min_value": 0.0,
            "negative_fraction": 0.0,
            "classification": "none",
        }
    
    min_val = float(np.min(integrand_real))
    negative_fraction = float(np.sum(negative_mask)) / integrand_real.size
    relative_magnitude = abs(min_val) / max(reference_scale, 1e-30)
    
    # Classify based on magnitude and spatial pattern
    if relative_magnitude < 1e-10:
        classification = "numerical_noise"
        severity = "info"
    elif relative_magnitude < 1e-6:
        classification = "physical_interference_expected"
        severity = "info"
    elif relative_magnitude < 1e-3:
        classification = "physical_interference_strong"
        severity = "warning" if ks < 0.3 else "info"
    else:
        classification = "calculation_error_suspected"
        severity = "error"
    
    # Check spatial coherence (physical interference should be coherent)
    if negative_fraction > 0.01:  # More than 1% negative points
        try:
            from scipy.ndimage import correlate
            # Compute autocorrelation of sign pattern
            sign_pattern = (integrand_real < 0).astype(float)
            autocorr = correlate(sign_pattern, sign_pattern, mode='constant')
            center = tuple(s // 2 for s in autocorr.shape)
            spatial_coherence = float(autocorr[center]) / sign_pattern.size
            
            if spatial_coherence < 0.3:
                classification += "_random_pattern_suspicious"
                severity = "warning"
        except ImportError:
            pass  # Skip spatial coherence check if scipy not available
    
    result = {
        "has_negatives": True,
        "min_value": min_val,
        "negative_fraction": negative_fraction,
        "relative_magnitude": relative_magnitude,
        "classification": classification,
    }
    
    if severity in ["warning", "error"]:
        message = (
            f"{label}: Negative integrand values detected.\n"
            f"  Min value: {min_val:.3e}\n"
            f"  Relative to peak: {relative_magnitude:.3e}\n"
            f"  Negative fraction: {negative_fraction:.2%}\n"
            f"  Classification: {classification}\n"
        )
        if classification == "calculation_error_suspected":
            message += (
                "  CRITICAL: Large negative values suggest calculation error.\n"
                "  Possible causes:\n"
                "    - Taking abs() before integration (destroys interference)\n"
                "    - Wrong sign in series coefficients\n"
                "    - Propagators evaluated at wrong spectral points\n"
            )
        elif "random_pattern" in classification:
            message += (
                "  WARNING: Random spatial pattern suggests numerical instability.\n"
                "  Physical interference should show coherent fringe patterns.\n"
            )
        _guard(False, message, severity=severity)
    
    return result


def validate_integration_domain(
    U: np.ndarray,
    V: np.ndarray,
    k: float,
    correlation_length: float,
    ks: float,
    label: str = "integration_domain",
) -> None:
    """Validate that integration domain is adequate for surface roughness.
    
    Based on re-derivation:
    - Main spectral peak: ~1/correlation_length
    - Interference side lobes: extend to ~5-10/correlation_length
    - Radiation condition: |u|, |v| < k for propagating modes
    
    Parameters
    ----------
    U, V : np.ndarray
        Spectral wavenumber grids
    k : float
        Incident wavenumber
    correlation_length : float
        Surface correlation length (meters)
    ks : float
        Normalized roughness kσ
    label : str
        Identifier
    """
    kl = k * correlation_length
    u_max = float(np.max(np.abs(U)))
    v_max = float(np.max(np.abs(V)))
    
    # Rule 1: Domain should extend to at least 3/kl (captures main lobe)
    min_extent = 3.0 / max(kl, 1e-6)
    _guard(
        u_max >= min_extent and v_max >= min_extent,
        f"{label}: Integration domain too small. "
        f"Max |u| = {u_max:.3f}, |v| = {v_max:.3f}, need ≥ {min_extent:.3f} (3/kl). "
        "Main spectral peak will be truncated.",
    )
    
    # Rule 2: For rough surfaces, should extend to 5-10/kl (captures side lobes)
    if ks > 0.5:
        recommended_extent = 7.0 / max(kl, 1e-6)  # Middle of range
        _guard(
            u_max >= recommended_extent and v_max >= recommended_extent,
            f"{label}: For kσ = {ks:.2f}, recommend domain ≥ {recommended_extent:.3f} (7/kl) "
            f"to capture interference side lobes. Current: {u_max:.3f}",
            severity="warning",
        )
    
    # Rule 3: Domain should not greatly exceed k (radiation condition)
    excessive_extent = 2.0 * k
    _guard(
        u_max <= excessive_extent or v_max <= excessive_extent,
        f"{label}: Integration domain excessively large. "
        f"Max |u| = {u_max:.3f}, |v| = {v_max:.3f} vs k = {k:.3f}. "
        "This wastes computation on evanescent modes.",
        severity="warning",
    )
    
    # Rule 4: Grid resolution should be fine enough
    du = float(U[1, 0] - U[0, 0]) if U.shape[0] > 1 else 0
    dv = float(V[0, 1] - V[0, 0]) if V.shape[1] > 1 else 0
    
    max_resolution = 0.5 / max(kl, 1e-6)  # Nyquist criterion
    _guard(
        abs(du) <= max_resolution and abs(dv) <= max_resolution,
        f"{label}: Grid too coarse. "
        f"du = {du:.4f}, dv = {dv:.4f}, recommend ≤ {max_resolution:.4f} (0.5/kl)",
        severity="warning",
    )


# =============================================================================
# SECTION 2: ENHANCED PHASE FACTOR VALIDATION
# =============================================================================


def validate_phase_factor_structure(
    phase_factor: complex,
    sigma: float,
    kz: float,
    ksz: float,
    q: float,
    qp: float,
    term_type: str,
    label: str,
) -> None:
    """Validate phase factor structure against Gaussian statistics theory.
    
    From re-derivation: Phase factors have form exp[-σ²·Ψ] where Ψ comes from
    moment generating function of Gaussian surface heights:
    
    For Kirchhoff-complementary (gkc1-3):
        Ψ = k²ₛz + k²z + kₛz·kz + q² - q·kₛz + q·kz
    
    For complementary (gc1-14):
        Ψ = k²ₛz + k²z + q² + q'² - (kₛz - kz)(q + q')
    
    These must satisfy:
    1. Ψ > 0 (ensures damping, not growth)
    2. Magnitude O(k²) (correct scaling)
    3. Smooth surface limit: exp[-σ²·Ψ] → 1 as σ → 0
    
    Parameters
    ----------
    phase_factor : complex
        The computed exp[-σ²·Ψ] value
    sigma : float
        RMS height
    kz, ksz : float
        Vertical wavenumbers
    q, qp : float
        Spectral vertical wavenumbers
    term_type : str
        'gkc' or 'gc' to select validation equations
    label : str
        Identifier
    """
    mag = abs(phase_factor)
    
    # Check 1: Magnitude should be ≤ 1 (damping, not growth)
    _guard(
        mag <= 1.0 + 1e-6,
        f"{label}: Phase factor |exp[-σ²·Ψ]| = {mag:.6f} > 1. "
        "This indicates exponential GROWTH instead of damping. "
        "Check sign in phase factor computation.",
    )
    
    # Check 2: Should not be exactly 1 unless σ is very small
    if sigma > 0.01:  # Non-negligible roughness
        _guard(
            mag < 0.999,
            f"{label}: Phase factor ≈ 1 despite σ = {sigma:.3f}. "
            "Height variance may not be included in phase calculation.",
            severity="warning",
        )
    
    # Check 3: Compute expected Ψ and verify
    if term_type == 'gkc':
        # Kirchhoff-complementary: Eq A1-A3 structure
        psi_expected = kz**2 + ksz**2 + kz*ksz + q**2 - q*ksz + q*kz
    else:
        # Complementary: Eq A4-A17 structure
        psi_expected = ksz**2 + kz**2 + q**2 + qp**2 - (ksz - kz)*(q + qp)
    
    # Check Ψ is positive
    _guard(
        psi_expected > -1e-9,
        f"{label}: Ψ = {psi_expected:.6f} ≤ 0. "
        f"Phase factor will cause exponential growth. "
        f"Check wavenumber signs: kz={kz:.3f}, ksz={ksz:.3f}, q={q:.3f}, qp={qp:.3f}",
    )
    
    # Check computed vs expected magnitude
    mag_expected = math.exp(-sigma**2 * psi_expected)
    relative_error = abs(mag - mag_expected) / max(mag_expected, 1e-30)
    
    _guard(
        relative_error < 0.1,
        f"{label}: Phase factor magnitude mismatch. "
        f"Computed: {mag:.6f}, Expected: {mag_expected:.6f} (error: {relative_error:.2%}). "
        "Check if Ψ formula matches paper equations.",
        severity="warning",
    )
    
    # Check 4: Imaginary part should be negligible (real exponential damping)
    imag_part = abs(np.imag(phase_factor))
    _guard(
        imag_part < 0.1 * mag,
        f"{label}: Phase factor has large imaginary part {imag_part:.3e}. "
        "Should be real-valued damping factor. Check for complex Ψ.",
        severity="warning",
    )


# =============================================================================
# SECTION 3: ENHANCED PROPAGATOR VALIDATION
# =============================================================================


def validate_propagator_shifts(
    U: np.ndarray,
    V: np.ndarray,
    U_shifted: np.ndarray,
    V_shifted: np.ndarray,
    kx: float,
    ky: float,
    ksx: float,
    ksy: float,
    term_index: int,
    label: str,
    *,
    tolerance: float = 1e-9,
) -> None:
    """Validate that propagators are evaluated at correct spectral points.
    """
    incident_shift_terms = {6, 7, 8, 9, 10, 11}
    scattered_shift_terms = {3, 4, 5, 12, 13, 14}
    if term_index == 1:
        expected_u = U
        expected_v = V
    elif term_index == 2:
        expected_u = -kx - ksx - U
        expected_v = -ky - ksy - V
    elif term_index in scattered_shift_terms:
        expected_u = -ksx - U
        expected_v = -ksy - V
    elif term_index in incident_shift_terms:
        expected_u = -kx - U
        expected_v = -ky - V
    else:
        # Not a recognized gc term; skip validation.
        return

    cond_u = np.allclose(U_shifted, expected_u, atol=tolerance, rtol=0.0)
    cond_v = np.allclose(V_shifted, expected_v, atol=tolerance, rtol=0.0)
    _guard(
        cond_u and cond_v,
        f"{label}: Propagator evaluated at incorrect coordinates for gc{term_index}.",
        severity="warning",
    )

    if term_index != 1 and U.size > 1 and V.size > 1:
        # Ensure the shifted coordinates vary across the mesh (not constant scalars)
        u_variation = float(np.std(U_shifted))
        v_variation = float(np.std(V_shifted))
        _guard(
            (u_variation > tolerance) or (v_variation > tolerance),
            f"{label}: Expected shifted coordinates to vary with integration point "
            f"but detected nearly constant values (σ_u={u_variation:.3e}, "
            f"σ_v={v_variation:.3e}).",
            severity="warning",
        )


def validate_propagator_conjugation(
    F_plus: np.ndarray,
    F_minus: np.ndarray,
    G_plus: np.ndarray,
    G_minus: np.ndarray,
    pol: str,
    label: str,
) -> None:
    """Validate proper conjugation handling in propagator products.
    
    From re-derivation: Equation 13 has structure:
        F⁺(u,v) F⁺*(u',v') gc_i(...)
    
    The * denotes complex conjugate. When computing:
        integrand = F⁺ × F⁺*
    
    We must use: F_plus * np.conj(F_plus), NOT F_plus * F_plus
    
    For same-point evaluation (i=1): |F⁺|² = F⁺ · F⁺* is real
    For different-point evaluation (i=2): F⁺(u,v) · F⁻*(u',v') is complex
    
    Parameters
    ----------
    F_plus, F_minus, G_plus, G_minus : np.ndarray
        Propagator arrays
    pol : str
        Polarization
    label : str
        Identifier
    """
    # Check that propagators are complex (or can be treated as such)
    _guard(
        np.iscomplexobj(F_plus) or np.all(np.imag(F_plus) == 0),
        f"{label}: F_plus should be complex or real-valued, got dtype {F_plus.dtype}",
    )
    
    # Check a sample product for proper conjugation
    # For same-point: |F⁺|² should be real and positive
    sample_product = F_plus * np.conj(F_plus)
    sample_real = np.real(sample_product)
    sample_imag = np.imag(sample_product)
    
    max_imag_ratio = np.max(np.abs(sample_imag)) / (np.max(np.abs(sample_real)) + 1e-30)
    _guard(
        max_imag_ratio < 0.01,
        f"{label}: |F⁺|² has significant imaginary part (ratio {max_imag_ratio:.3e}). "
        "Product F × F* should be real. Check conjugation implementation.",
        severity="warning",
    )
    
    # Check that products are not obviously wrong (e.g., using F*F instead of F*F*)
    wrong_product = F_plus * F_plus  # Wrong: no conjugate
    right_product = F_plus * np.conj(F_plus)  # Right: with conjugate
    
    # They should differ significantly if F_plus is complex
    if np.iscomplexobj(F_plus):
        relative_diff = (
            np.mean(np.abs(wrong_product - right_product)) / 
            (np.mean(np.abs(right_product)) + 1e-30)
        )
        _guard(
            relative_diff > 1e-6,
            f"{label}: F×F and F×F* are nearly identical. "
            "This suggests propagators are purely real (missing physics) or "
            "conjugation is not being applied.",
            severity="warning",
        )


def validate_spectral_coupling_structure(
    spectrum_1: np.ndarray,
    spectrum_2: np.ndarray,
    k1_x: float,
    k1_y: float,
    k2_x: float,
    k2_y: float,
    label: str,
) -> None:
    """Validate spectral coupling W^(m)(k₁) × W^(n)(k₂) structure.
    
    From re-derivation: Multiple scattering involves products of spectra
    evaluated at DIFFERENT spectral points:
    
    - W^(m)(kx + u, ky + v): First scattering event
    - W^(n)(ksx + u, ksy + v): Second scattering event
    
    These should:
    1. Both be non-negative (power spectrum property)
    2. Peak at different locations (unless backscatter)
    3. Have correlation if events are spatially close
    
    Parameters
    ----------
    spectrum_1, spectrum_2 : np.ndarray
        The two spectrum arrays
    k1_x, k1_y, k2_x, k2_y : float
        Center wavenumbers for each spectrum
    label : str
        Identifier
    """
    # Check 1: Both should be non-negative
    min_1 = np.min(spectrum_1)
    min_2 = np.min(spectrum_2)
    
    _guard(
        min_1 >= -1e-10,
        f"{label}: Spectrum 1 has negative values (min = {min_1:.3e}). "
        "Power spectrum must be non-negative.",
    )
    
    _guard(
        min_2 >= -1e-10,
        f"{label}: Spectrum 2 has negative values (min = {min_2:.3e}). "
        "Power spectrum must be non-negative.",
    )
    
    # Check 2: Product should preserve non-negativity
    product = spectrum_1 * spectrum_2
    min_product = np.min(product)
    
    _guard(
        min_product >= -1e-10,
        f"{label}: Spectrum product has negative values (min = {min_product:.3e}). "
        "Check if spectra are computed correctly.",
    )
    
    # Check 3: For backscatter, k1 and k2 should be related
    # k1 = kx + u, k2 = ksx + u
    # For backscatter: ksx = -kx, so k2 = -kx + u = -(kx - u)
    if abs(k1_x + k2_x) < 0.1 and abs(k1_y + k2_y) < 0.1:
        # Backscatter geometry detected
        # Spectra should be related by symmetry
        # (This is a heuristic check, not strict)
        pass
    
    # Check 4: Both spectra should have significant overlap in spectral domain
    # (Otherwise integration will be inefficient)
    overlap = np.sum((spectrum_1 > 1e-10) & (spectrum_2 > 1e-10))
    total_points = spectrum_1.size
    overlap_fraction = overlap / total_points
    
    _guard(
        overlap_fraction > 0.05,
        f"{label}: Spectra have < 5% overlap ({overlap_fraction:.1%}). "
        "This suggests spectral arguments may be incorrect, causing "
        "inefficient integration.",
        severity="warning",
    )


# =============================================================================
# SECTION 4: ENHANCED SERIES CONVERGENCE VALIDATION
# =============================================================================


def validate_double_summation_structure(
    n_max: int,
    m_max: int,
    series_terms_n: np.ndarray,
    series_terms_m: np.ndarray,
    label: str,
) -> None:
    """Validate double summation Σ_m Σ_n structure for multiple scattering.
    
    From re-derivation: Multiple scattering involves products:
        Σ_m [a_m]^m/m! W^(m)(k₁) × Σ_n [a_n]^n/n! W^(n)(k₂)
    
    Both series must converge independently. Check:
    1. Terms decrease monotonically (mostly)
    2. Ratio |term_{n+1}/term_n| < 1 for large n
    3. Final terms are small relative to sum
    
    Parameters
    ----------
    n_max, m_max : int
        Maximum summation indices
    series_terms_n, series_terms_m : np.ndarray
        Array of term magnitudes [term_1, term_2, ..., term_nmax]
    label : str
        Identifier
    """
    # Validate n-series
    if len(series_terms_n) > 2:
        # Check term ratios
        ratios_n = np.abs(series_terms_n[1:] / np.maximum(series_terms_n[:-1], 1e-30))
        
        # For convergent series, ratios should eventually be < 1
        late_ratios_n = ratios_n[-min(5, len(ratios_n)):]
        _guard(
            np.median(late_ratios_n) < 1.0,
            f"{label}: n-series not converging. "
            f"Late term ratios: {late_ratios_n}. "
            f"Expected median < 1.0 for convergent series.",
        )
        
        # Check final term is small
        if len(series_terms_n) > 0:
            final_ratio_n = series_terms_n[-1] / (np.sum(series_terms_n) + 1e-30)
            _guard(
                final_ratio_n < 0.01,
                f"{label}: n-series final term ({series_terms_n[-1]:.3e}) is "
                f"{final_ratio_n:.1%} of total sum. Recommend n_max ≥ {n_max + 10}",
                severity="warning",
            )
    
    # Validate m-series (same checks)
    if len(series_terms_m) > 2:
        ratios_m = np.abs(series_terms_m[1:] / np.maximum(series_terms_m[:-1], 1e-30))
        late_ratios_m = ratios_m[-min(5, len(ratios_m)):]
        
        _guard(
            np.median(late_ratios_m) < 1.0,
            f"{label}: m-series not converging. "
            f"Late term ratios: {late_ratios_m}.",
        )
        
        if len(series_terms_m) > 0:
            final_ratio_m = series_terms_m[-1] / (np.sum(series_terms_m) + 1e-30)
            _guard(
                final_ratio_m < 0.01,
                f"{label}: m-series final term ({series_terms_m[-1]:.3e}) is "
                f"{final_ratio_m:.1%} of total. Recommend m_max ≥ {m_max + 10}",
                severity="warning",
            )


def validate_series_term_ratios(
    terms: np.ndarray,
    coefficients: np.ndarray,
    sigma: float,
    label: str,
) -> None:
    """Validate series term ratios against theoretical expectations.
    
    From re-derivation: Series has form Σ [σ²·A]^n/n! where A ~ O(k²)
    
    Term ratio: |term_{n+1}/term_n| = |σ²·A/(n+1)|
    
    For convergence: |σ²·A| < n for large n
    For exponential convergence: |σ²·A| << n
    
    Parameters
    ----------
    terms : np.ndarray
        Computed term magnitudes
    coefficients : np.ndarray
        The 'A' factors for each term
    sigma : float
        RMS height
    label : str
        Identifier
    """
    if len(terms) < 3:
        return
    
    # Compute actual ratios
    actual_ratios = np.abs(terms[1:] / np.maximum(terms[:-1], 1e-30))
    
    # Compute theoretical ratios: σ²|A|/(n+1)
    theoretical_ratios = []
    for n in range(len(coefficients) - 1):
        A_magnitude = abs(coefficients[n])
        expected_ratio = (sigma**2 * A_magnitude) / (n + 2)  # n+2 because n is 0-indexed
        theoretical_ratios.append(expected_ratio)
    
    theoretical_ratios = np.array(theoretical_ratios[:len(actual_ratios)])
    
    # Compare
    if len(theoretical_ratios) > 0:
        # They should agree within factor of ~2 (due to W^(n) spectrum variations)
        ratio_agreement = actual_ratios / (theoretical_ratios + 1e-30)
        median_agreement = np.median(ratio_agreement)
        
        _guard(
            0.1 < median_agreement < 10.0,
            f"{label}: Series term ratios disagree with theory by factor {median_agreement:.2f}. "
            "Check if σ² factor is correctly included in series coefficients.",
            severity="warning",
        )


def detect_series_divergence(
    terms: np.ndarray,
    label: str,
    threshold_ratio: float = 1.5,
) -> None:
    """Detect if series is diverging (terms increasing).
    
    Divergent series indicates:
    1. Surface too rough for perturbation theory (kσ > 3)
    2. Wrong sign in series coefficients
    3. Missing damping factors
    
    Parameters
    ----------
    terms : np.ndarray
        Term magnitudes
    label : str
        Identifier
    threshold_ratio : float
        Ratio above which we consider series diverging
    """
    if len(terms) < 5:
        return
    
    # Check last few terms for growth
    late_terms = terms[-5:]
    if np.all(late_terms > 1e-30):
        growth_ratios = late_terms[1:] / late_terms[:-1]
        
        if np.any(growth_ratios > threshold_ratio):
            max_growth = np.max(growth_ratios)
            _guard(
                False,
                f"{label}: SERIES DIVERGENCE detected! "
                f"Terms increasing with ratio {max_growth:.2f} (threshold {threshold_ratio}). "
                "Possible causes:\n"
                "  1. Surface too rough for perturbation theory\n"
                "  2. Wrong sign in series coefficients\n"
                "  3. Missing exponential damping factor\n"
                "  4. Numerical instability",
                severity="error",
            )


# =============================================================================
# SECTION 5: BASIC INPUT VALIDATION (from original)
# =============================================================================


def validate_inputs(
    eps_r: complex,
    wavelength: float,
    sigma_m: float,
    corr_length_m: float,
    theta_i: float,
    theta_s: float,
) -> None:
    """Validate basic physical inputs."""
    _guard(wavelength > 0.0, "Vacuum wavelength must be positive")
    _guard(sigma_m > 0.0, "Surface RMS height must be positive")
    _guard(corr_length_m > 0.0, "Correlation length must be positive")
    _guard(
        eps_r.real >= 1.0 - 1e-9,
        f"Relative permittivity real part must be ≥ 1, got {eps_r.real!r}.",
    )
    _guard(
        eps_r.imag >= -1e-9,
        f"Relative permittivity imaginary part must be ≥ 0, got {eps_r.imag!r}.",
    )
    _guard(
        0.0 <= theta_i < 0.5 * math.pi,
        f"Incident angle must be in [0, π/2), got {theta_i!r}",
    )
    _guard(
        0.0 <= theta_s < 0.5 * math.pi,
        f"Scattered angle must be in [0, π/2), got {theta_s!r}",
    )


def validate_wavenumbers(k: float, theta_i: float, theta_s: float) -> None:
    """Ensure vertical wavenumber branches satisfy the radiation condition."""
    kiz = k * math.cos(theta_i)
    ksz = k * math.cos(theta_s)
    
    _guard(
        np.real(kiz) >= -1e-9,
        f"Incident vertical wavenumber k_z must be ≥ 0, got {kiz!r}",
    )
    _guard(
        np.real(ksz) >= -1e-9,
        f"Scattered vertical wavenumber k_sz must be ≥ 0, got {ksz!r}",
    )


def validate_radiation_condition(
    q1: complex,
    q2: complex,
    eps_r: complex,
    label: str = "spectral_point",
) -> None:
    """Validate radiation condition for spectral integration points."""
    # Upper medium: should be real for propagating modes
    if np.imag(q1) > 1e-10:
        warnings.warn(
            f"{label}: Upper medium q1 has imaginary part {np.imag(q1):.3e}. "
            "This is an evanescent mode (non-propagating).",
            UserWarning,
        )
    
    # Lower medium: MUST have positive imaginary part for lossy media
    _guard(
        np.imag(q2) >= -1e-12,
        f"{label}: Lower medium q2 has Im(q2) = {np.imag(q2):.3e} < 0. "
        "Check square root branch cut in q2 = √(εr·k² - u² - v²).",
    )


def validate_fresnel_bounds(Rv: complex, Rh: complex, eps_r: complex) -> None:
    """Check Fresnel reflection magnitudes stay within physical bounds."""
    max_mag = max(abs(Rv), abs(Rh))
    _guard(
        max_mag <= 1.0 + 1e-9,
        f"Fresnel magnitude exceeds unity: max(|Rv|, |Rh|) = {max_mag!r}. "
        "Check permittivity or angle calculation.",
    )


def validate_field_coefficients(coeffs: Mapping[str, complex]) -> None:
    """Check Kirchhoff field coefficients are finite and moderate in magnitude."""
    for label, value in coeffs.items():
        _guard(
            np.isfinite(value),
            f"Kirchhoff coefficient {label} is not finite: {value!r}",
        )
        _guard(
            abs(value) < 100.0,
            f"Kirchhoff coefficient {label} magnitude {abs(value):.3e} is "
            "unreasonably large (expect O(1)).",
            severity="warning",
        )


def ensure_array_finite(label: str, array: np.ndarray) -> None:
    """Ensure numerical arrays do not contain NaN or Inf."""
    _guard(
        np.all(np.isfinite(array)),
        f"Numerical array '{label}' contains non-finite values (NaN or Inf).",
    )


def validate_sigma_real(label: str, value) -> None:
    """Ensure scattering coefficients are (numerically) real and non-negative."""
    arr = np.asarray(value)
    if np.iscomplexobj(arr):
        real_val = float(np.real(arr))
        imag_val = float(np.imag(arr))
        scale = max(1.0, abs(real_val))
        imag_bound = 1e-12 * scale + 1e-15
        _guard(
            abs(imag_val) <= imag_bound,
            f"{label} scattering coefficient has large imaginary part: {value!r}.",
        )
        real_part = real_val
    else:
        real_part = float(arr)

    # Allow small negative values due to numerical precision
    tolerance = 1e-6
    _guard(
        real_part >= -tolerance,
        f"{label} scattering coefficient significantly negative: {real_part:.3e}.",
    )


# =============================================================================
# SECTION 6: MULTIPLE SCATTERING SPECIFIC VALIDATION
# =============================================================================


def validate_ms_balance(
    sigma_single: float,
    sigma_multiple: float,
    ks: float,
    polarization: str,
) -> None:
    """Ensure multiple scattering stays bounded relative to single scattering."""
    pol = polarization.upper()
    is_copol = pol in {"HH", "VV"}
    
    if not is_copol:
        return
    
    # CRITICAL: Multiple should NEVER exceed single for co-pol
    if sigma_multiple > sigma_single * 1.01:
        ratio = sigma_multiple / sigma_single
        ratio_db = 10.0 * math.log10(ratio)
        _guard(
            False,
            f"*** CRITICAL VIOLATION *** Co-pol {pol}: Multiple scattering "
            f"({sigma_multiple:.3e}) EXCEEDS single scattering ({sigma_single:.3e}). "
            f"Ratio = {ratio:.3f} ({ratio_db:+.2f} dB). kσ = {ks:.3f}.",
            severity="error",
        )


def validate_cross_pol_single_scattering(
    sigma_s_hv: float,
    sigma_m_hv: float,
    *,
    threshold_ratio: float = 0.01,
) -> None:
    """Cross-pol comes ONLY from multiple scattering."""
    if sigma_m_hv <= 0.0:
        return
    
    if sigma_s_hv <= 0.0:
        return
    
    ratio = sigma_s_hv / sigma_m_hv
    
    _guard(
        ratio < threshold_ratio,
        f"*** CRITICAL VIOLATION *** Cross-pol single scattering non-zero: "
        f"σ°_hv(single) = {sigma_s_hv:.3e}, σ°_hv(multiple) = {sigma_m_hv:.3e}, "
        f"ratio = {ratio:.3e} > {threshold_ratio}.",
    )


def validate_kirchhoff_polarization_independence(
    f_hh: complex,
    f_vv: complex,
    *,
    rtol: float = 0.20,
) -> None:
    """Verify Kirchhoff coefficient is polarization-independent."""
    if abs(f_hh) < 1e-12 and abs(f_vv) < 1e-12:
        return
    
    diff = abs(f_hh - f_vv)
    ref = max(abs(f_hh), abs(f_vv))
    
    _guard(
        diff <= rtol * ref,
        f"Kirchhoff polarization independence: "
        f"f_hh = {f_hh:.6e}, f_vv = {f_vv:.6e}, diff = {diff:.6e}.",
        severity="warning",
    )


def validate_cross_pol_kirchhoff_zero(
    sigma_k_hv: float,
    *,
    threshold_db: float = -80.0,
) -> None:
    """Verify Kirchhoff cross-pol vanishes in backscatter."""
    if sigma_k_hv <= 0.0:
        return
    
    sigma_k_hv_db = 10.0 * math.log10(max(sigma_k_hv, 1e-30))
    
    _guard(
        sigma_k_hv_db < threshold_db,
        f"Kirchhoff cross-pol in backscatter: "
        f"σ°^k_hv = {sigma_k_hv:.3e} ({sigma_k_hv_db:.1f} dB) > {threshold_db} dB.",
    )


def validate_component_hierarchy(
    sigma_k: float,
    sigma_kc: float,
    sigma_c: float,
    ks: float,
    pol: str,
) -> None:
    """Validate single scattering component magnitude hierarchy."""
    is_copol = pol.upper() in {"HH", "VV"}
    
    if not is_copol:
        return
    
    if ks < 0.3:
        _guard(
            sigma_k > sigma_kc * 0.1,
            f"Co-pol {pol}: Kirchhoff ({sigma_k:.3e}) should dominate cross term "
            f"({sigma_kc:.3e}) for kσ = {ks:.3f} < 0.3.",
            severity="warning",
        )


def validate_hh_vv_ms_asymmetry(
    sigma_s_hh: float,
    sigma_m_hh: float,
    sigma_s_vv: float,
    sigma_m_vv: float,
    ks: float,
) -> None:
    """Validate HH vs VV multiple scattering asymmetry."""
    if sigma_s_hh <= 0.0 or sigma_s_vv <= 0.0:
        return
    
    if sigma_m_hh <= 0.0 or sigma_m_vv <= 0.0:
        return
    
    ms_impact_hh = sigma_m_hh / sigma_s_hh
    ms_impact_vv = sigma_m_vv / sigma_s_vv
    
    if ks > 0.3:
        _guard(
            ms_impact_hh >= ms_impact_vv * 0.8,
            f"Multiple scattering impact: HH {ms_impact_hh:.3f}, "
            f"VV {ms_impact_vv:.3f} for kσ = {ks:.3f}.",
            severity="warning",
        )


# =============================================================================
# SECTION 7: PHYSICAL BOUNDS AND CONSERVATION LAWS
# =============================================================================


def ensure_reciprocity(hv: float, vh: float, *, rtol: float = 1e-6) -> None:
    """Check monostatic reciprocity σ_hv == σ_vh."""
    if hv == 0.0 and vh == 0.0:
        return
    diff = abs(hv - vh)
    ref = max(abs(hv), abs(vh), 1e-12)
    _guard(
        diff <= rtol * ref,
        f"Reciprocity violated: σ_hv = {hv:.6e}, σ_vh = {vh:.6e}",
    )


def validate_energy_conservation(
    sigma_hh: float,
    sigma_vv: float,
    sigma_hv: float,
    *,
    threshold: float = 1.0,
    severity: str = "error",
) -> None:
    """Validate energy conservation: all scattering coefficients ≤ 1."""
    for pol, sigma in [("HH", sigma_hh), ("VV", sigma_vv), ("HV", sigma_hv)]:
        _guard(
            sigma <= threshold + 1e-6,
            f"*** ENERGY CONSERVATION VIOLATED *** σ°_{pol} = {sigma:.3e} > 1.",
            severity=severity,
        )


def validate_cross_pol_ordering(
    sigma_hv: float,
    sigma_hh: float,
    sigma_vv: float,
) -> None:
    """Validate cross-pol magnitude ordering."""
    min_copol = min(sigma_hh, sigma_vv)

    if sigma_hv > min_copol * 2.0:
        warnings.warn(
            f"Cross-pol σ°_hv = {sigma_hv:.3e} exceeds 2× minimum co-pol "
            f"({min_copol:.3e}).",
            UserWarning,
        )


def validate_cross_pol_magnitude(
    sigma_hv: float,
    ks: float,
    *,
    max_reasonable_db: float = 10.0,
) -> None:
    """Validate cross-pol has reasonable magnitude.

    Cross-pol backscatter should be significantly weaker than unity
    (0 dB) for natural surfaces. Values above ~10 dB (linear = 10)
    suggest computational errors.
    """
    if sigma_hv <= 0.0:
        return

    sigma_hv_db = 10.0 * math.log10(sigma_hv)

    _guard(
        sigma_hv_db < max_reasonable_db,
        f"*** CRITICAL VIOLATION *** Cross-pol magnitude unreasonably large: "
        f"σ°_hv = {sigma_hv:.3e} ({sigma_hv_db:.1f} dB) > {max_reasonable_db} dB. "
        f"kσ = {ks:.3f}. This likely indicates a computational error.",
        severity="error",
    )


def validate_nadir_symmetry(
    sigma_hh: float,
    sigma_vv: float,
    theta: float,
    *,
    rtol: float = 0.05,
) -> None:
    """Validate HH = VV at normal incidence."""
    if theta > 5.0 * math.pi / 180.0:
        return
    
    if sigma_hh <= 0.0 or sigma_vv <= 0.0:
        return
    
    ratio = abs(sigma_hh - sigma_vv) / max(sigma_hh, sigma_vv)
    
    _guard(
        ratio <= rtol,
        f"Nadir symmetry violation at θ = {math.degrees(theta):.1f}°: "
        f"σ°_hh = {sigma_hh:.3e}, σ°_vv = {sigma_vv:.3e}",
        severity="warning",
    )


def validate_brewster_behavior(
    sigma_vv: float,
    sigma_hh: float,
    theta: float,
    eps_r: complex,
    *,
    tolerance_deg: float = 5.0,
) -> None:
    """Check for Brewster angle minimum in VV (not HH)."""
    theta_deg = math.degrees(theta)
    brewster_deg = math.degrees(math.atan(math.sqrt(eps_r.real)))
    
    if abs(theta_deg - brewster_deg) > tolerance_deg:
        return
    
    if sigma_vv > 0.0 and sigma_hh > 0.0:
        ratio = sigma_vv / sigma_hh
        _guard(
            ratio < 0.9,
            f"Near Brewster angle (θ = {theta_deg:.1f}°, θ_B = {brewster_deg:.1f}°): "
            f"σ°_vv/σ°_hh = {ratio:.3f} ≥ 0.9.",
            severity="warning",
        )


def validate_angular_trends(
    sigma_at_angles: list[Tuple[float, float]],
    pol: str,
    *,
    check_monotonic: bool = True,
) -> None:
    """Validate angular behavior trends."""
    if len(sigma_at_angles) < 3:
        return
    
    angles = [theta for theta, _ in sigma_at_angles]
    sigmas = [sigma for _, sigma in sigma_at_angles]
    
    if check_monotonic and pol.upper() in {"HH"}:
        increases = 0
        for i in range(1, len(sigmas)):
            if sigmas[i] > sigmas[i-1] * 1.1:
                increases += 1
        
        if increases > len(sigmas) // 3:
            warnings.warn(
                f"{pol}: Scattering coefficient increases with angle at "
                f"{increases}/{len(sigmas)} points.",
                UserWarning,
            )


# =============================================================================
# SECTION 8: INTEGRATION AND CONVERGENCE
# =============================================================================


def assert_spectrum_consistency(k: float, corr_length_m: float, kl: float) -> None:
    """Verify that kl and correlation length use consistent units."""
    expected = k * corr_length_m
    _guard(
        math.isclose(expected, kl, rel_tol=1e-4, abs_tol=1e-6),
        f"Inconsistent correlation length units: expected kl = {expected:.6e}, got {kl:.6e}.",
    )


def enforce_nonnegative_integrand(
    label: str,
    integrand: np.ndarray,
    *,
    tolerance: float = 1e-5,
) -> None:
    """Ensure complementary integrands remain non-negative within tolerance."""
    min_val = float(np.min(np.real(integrand)))
    _guard(
        min_val >= -tolerance,
        f"{label} integrand has negative values: min = {min_val:.3e} < -{tolerance}.",
    )


def enforce_min_integrand_energy(
    label: str,
    integrand: np.ndarray,
    reference_scale: float,
    *,
    min_ratio: float = 1e-6,
) -> None:
    """Ensure complementary integrands retain sufficient magnitude."""
    _guard(
        reference_scale > 0.0,
        f"{label} reference scale must be positive, got {reference_scale}",
    )
    
    mean_mag = float(np.mean(np.abs(integrand)))
    threshold = min_ratio * reference_scale
    
    _guard(
        mean_mag >= threshold,
        f"{label} integrand energy collapsed: mean|I| = {mean_mag:.3e}, "
        f"expected ≥ {threshold:.3e}.",
    )


def validate_series_convergence(
    terms: np.ndarray,
    label: str,
    *,
    threshold: float = 1e-4,
) -> None:
    """Validate series convergence by checking term ratios."""
    if len(terms) < 3:
        return
    
    for i in range(len(terms) - 3, len(terms) - 1):
        if abs(terms[i]) < 1e-30:
            continue
        
        ratio = abs(terms[i+1] / terms[i])
        if ratio > 1.0:
            warnings.warn(
                f"{label}: Series term {i+1} larger than term {i} "
                f"(ratio {ratio:.3e} > 1).",
                UserWarning,
            )
    
    total = abs(np.sum(terms))
    if total > 0:
        final_ratio = abs(terms[-1]) / total
        _guard(
            final_ratio < threshold,
            f"{label}: Final term {abs(terms[-1]):.3e} is {final_ratio:.3e} "
            f"of total. Series not converged.",
            severity="warning",
        )


def validate_phase_factors(
    phase_factor: complex,
    label: str,
    ks: float,
) -> None:
    """Validate exponential phase/damping factors."""
    mag = abs(phase_factor)
    
    _guard(
        mag <= 1.0 + 1e-6,
        f"{label}: Phase factor magnitude {mag:.3e} > 1.",
    )
    
    if ks > 0.3 and mag > 0.99:
        warnings.warn(
            f"{label}: Phase factor magnitude {mag:.6f} ≈ 1 for kσ = {ks:.3f}.",
            UserWarning,
        )


def validate_spectral_coupling(
    has_product_spectra: bool,
    is_multiple_scattering: bool,
) -> None:
    """Verify spectral coupling structure for multiple scattering."""
    if is_multiple_scattering:
        _guard(
            has_product_spectra,
            "Multiple scattering term does not have product of spectra.",
        )


def validate_double_summation_convergence(
    n_terms_used: int,
    m_terms_used: int,
    ks: float,
    *,
    min_terms: int = 5,
) -> None:
    """Validate double summation convergence for multiple scattering."""
    _guard(
        n_terms_used >= min_terms,
        f"Double summation n_terms = {n_terms_used} < {min_terms}.",
    )
    
    _guard(
        m_terms_used >= min_terms,
        f"Double summation m_terms = {m_terms_used} < {min_terms}.",
    )
    
    if ks > 0.5:
        recommended = int(20 + 20 * (ks - 0.5))
        if n_terms_used < recommended or m_terms_used < recommended:
            warnings.warn(
                f"For kσ = {ks:.3f}, recommend n_max = m_max ≥ {recommended}.",
                UserWarning,
            )


# =============================================================================
# SECTION 9: FIELD COEFFICIENTS
# =============================================================================


def validate_field_coefficient_structure(
    F_coeffs: Mapping[str, complex],
    G_coeffs: Mapping[str, complex],
    eps_r: complex,
    Rv: complex,
    Rh: complex,
) -> None:
    """Validate complementary field coefficient structure."""
    for label, F in F_coeffs.items():
        _guard(
            np.isfinite(F),
            f"Complementary coefficient F_{label} is not finite: {F!r}",
        )
    
    for label, G in G_coeffs.items():
        _guard(
            np.isfinite(G),
            f"Complementary coefficient G_{label} is not finite: {G!r}",
        )
    
    if abs(eps_r) > 1000.0:
        max_G_mag = max(abs(G) for G in G_coeffs.values()) if G_coeffs else 0.0
        _guard(
            max_G_mag < 0.1,
            f"Near-PEC limit but max|G±| = {max_G_mag:.3e}.",
            severity="warning",
        )


def validate_complementary_coefficients(
    F_plus: complex,
    F_minus: complex,
    G_plus: complex,
    G_minus: complex,
    pol: str,
) -> None:
    """Check F±, G± satisfy expected relationships."""
    if abs(F_plus) > 1e-12 and abs(F_minus) > 1e-12:
        ratio = abs((F_plus - F_minus) / max(abs(F_plus), abs(F_minus)))
        _guard(
            ratio > 1e-6,
            f"Complementary {pol}: F+ and F- are nearly identical.",
            severity="warning",
        )


# =============================================================================
# SECTION 10: LIMITING CASES
# =============================================================================


def validate_pec_limit(
    sigma_computed: float,
    sigma_pec_analytical: float,
    eps_r: complex,
    pol: str,
    *,
    threshold_db: float = 0.5,
) -> None:
    """Validate PEC limiting case."""
    if abs(eps_r) < 100.0:
        return
    
    if sigma_pec_analytical <= 0.0:
        return
    
    diff_db = 10.0 * math.log10(abs(sigma_computed / sigma_pec_analytical))
    
    _guard(
        abs(diff_db) < threshold_db,
        f"PEC limit validation failed for {pol}: diff = {diff_db:.2f} dB.",
        severity="warning",
    )


def validate_spm_regime(
    sigma_aiem: float,
    sigma_spm: float,
    ks: float,
    rms_slope: float,
    pol: str,
    *,
    threshold_db: float = 2.0,
) -> None:
    """Validate Small Perturbation Method regime."""
    if ks >= 0.3 or rms_slope >= 0.3:
        return
    
    if sigma_spm <= 0.0:
        return
    
    diff_db = 10.0 * math.log10(abs(sigma_aiem / sigma_spm))
    
    _guard(
        abs(diff_db) < threshold_db,
        f"SPM regime validation failed for {pol}: diff = {diff_db:.2f} dB.",
        severity="warning",
    )


def validate_smooth_limit(
    sigma_total: float,
    sigma_specular_expected: float,
    ks: float,
    pol: str,
) -> None:
    """Validate smooth surface limit."""
    if ks > 0.1:
        return
    
    is_crosspol = pol.upper() in {"HV", "VH"}
    
    if is_crosspol:
        sigma_db = 10.0 * math.log10(max(sigma_total, 1e-30))
        _guard(
            sigma_db < -40.0,
            f"Smooth limit: Cross-pol {pol} = {sigma_db:.1f} dB > -40 dB.",
            severity="warning",
        )


# =============================================================================
# SECTION 11: COMPREHENSIVE VALIDATION SUITES
# =============================================================================


def validate_aiem_single_scattering(
    sigma_k_hh: float,
    sigma_k_vv: float,
    sigma_k_hv: float,
    sigma_kc_hh: float,
    sigma_kc_vv: float,
    sigma_kc_hv: float,
    sigma_c_hh: float,
    sigma_c_vv: float,
    sigma_c_hv: float,
    f_hh: complex,
    f_vv: complex,
    ks: float,
) -> None:
    """Comprehensive validation of single scattering components."""
    validate_kirchhoff_polarization_independence(f_hh, f_vv)
    validate_cross_pol_kirchhoff_zero(sigma_k_hv, threshold_db=-80.0)
    validate_component_hierarchy(sigma_k_hh, sigma_kc_hh, sigma_c_hh, ks, "HH")
    validate_component_hierarchy(sigma_k_vv, sigma_kc_vv, sigma_c_vv, ks, "VV")
    
    sigma_s_hh = sigma_k_hh + sigma_kc_hh + sigma_c_hh
    sigma_s_vv = sigma_k_vv + sigma_kc_vv + sigma_c_vv
    sigma_s_hv = sigma_k_hv + sigma_kc_hv + sigma_c_hv
    
    validate_energy_conservation(
        sigma_s_hh, sigma_s_vv, sigma_s_hv,
        threshold=1.0, severity="warning"
    )


def validate_aiem_multiple_scattering(
    sigma_m_hh: float,
    sigma_m_vv: float,
    sigma_m_hv: float,
    sigma_s_hh: float,
    sigma_s_vv: float,
    sigma_s_hv: float,
    ks: float,
    has_double_summation: bool,
    n_terms: int,
    m_terms: int,
) -> None:
    """Comprehensive validation of multiple scattering."""
    validate_cross_pol_single_scattering(sigma_s_hv, sigma_m_hv, threshold_ratio=0.01)
    validate_hh_vv_ms_asymmetry(sigma_s_hh, sigma_m_hh, sigma_s_vv, sigma_m_vv, ks)
    validate_ms_balance(sigma_s_hh, sigma_m_hh, ks, "HH")
    validate_ms_balance(sigma_s_vv, sigma_m_vv, ks, "VV")
    validate_spectral_coupling(has_double_summation, is_multiple_scattering=True)
    validate_double_summation_convergence(n_terms, m_terms, ks)
    validate_energy_conservation(
        sigma_m_hh, sigma_m_vv, sigma_m_hv,
        threshold=1.0, severity="warning"
    )


def validate_aiem_multiple_scattering_comprehensive(
    sigma_m_hh: float,
    sigma_m_vv: float,
    sigma_m_hv: float,
    sigma_s_hh: float,
    sigma_s_vv: float,
    sigma_s_hv: float,
    ks: float,
    n_max: int,
    m_max: int,
    integrand_kc_hh: np.ndarray,
    integrand_c_hh: np.ndarray,
    integrand_kc_hv: np.ndarray,
    integrand_c_hv: np.ndarray,
    U: np.ndarray,
    V: np.ndarray,
    k: float,
    correlation_length: float,
    propagators: Dict[str, np.ndarray],
) -> None:
    """Comprehensive validation suite for multiple scattering implementation."""
    # Integration domain
    validate_integration_domain(U, V, k, correlation_length, ks, "MS_domain")
    
    # Integrand structure for HH
    validate_integrand_structure(
        integrand_kc_hh + integrand_c_hh,
        "MS_HH_total",
        expected_peak_location=None,
        ks=ks,
    )
    
    # Negative value classification
    reference_hh = np.max(np.abs(integrand_kc_hh + integrand_c_hh))
    neg_class_hh = classify_negative_integrand(
        integrand_kc_hh + integrand_c_hh,
        "MS_HH",
        reference_hh,
        ks,
    )
    
    # HV integrand
    if sigma_m_hv > 0:
        reference_hv = np.max(np.abs(integrand_kc_hv + integrand_c_hv))
        neg_class_hv = classify_negative_integrand(
            integrand_kc_hv + integrand_c_hv,
            "MS_HV",
            reference_hv,
            ks,
        )
        
        if ks > 0.3 and not neg_class_hv["has_negatives"]:
            warnings.warn(
                "MS_HV: No negative values detected for kσ > 0.3.",
                UserWarning,
            )
    
    # MS balance
    validate_ms_balance(sigma_s_hh, sigma_m_hh, ks, "HH")
    validate_ms_balance(sigma_s_vv, sigma_m_vv, ks, "VV")
    
    # Energy conservation
    sigma_total_hh = sigma_s_hh + sigma_m_hh
    sigma_total_vv = sigma_s_vv + sigma_m_vv
    sigma_total_hv = sigma_s_hv + sigma_m_hv
    
    validate_energy_conservation(
        sigma_total_hh, sigma_total_vv, sigma_total_hv,
        threshold=1.0, severity="warning",
    )


def validate_aiem_total_scattering(
    sigma_hh: float,
    sigma_vv: float,
    sigma_hv: float,
    sigma_vh: float,
    sigma_s_hh: float,
    sigma_s_vv: float,
    sigma_m_hh: float,
    sigma_m_vv: float,
    theta: float,
    eps_r: complex,
) -> None:
    """Comprehensive validation of total scattering."""
    validate_energy_conservation(
        sigma_hh, sigma_vv, sigma_hv, threshold=1.0, severity="error"
    )
    ensure_reciprocity(sigma_hv, sigma_vh, rtol=1e-6)
    validate_cross_pol_ordering(sigma_hv, sigma_hh, sigma_vv)
    validate_nadir_symmetry(sigma_hh, sigma_vv, theta, rtol=0.05)
    validate_brewster_behavior(sigma_vv, sigma_hh, theta, eps_r)
    
    _guard(
        abs(sigma_hh - (sigma_s_hh + sigma_m_hh)) / max(sigma_hh, 1e-30) < 0.01,
        f"Total HH decomposition error",
    )
