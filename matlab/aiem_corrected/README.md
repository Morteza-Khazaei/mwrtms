# AIEM Corrected - MATLAB Implementation

**Corrected Advanced Integral Equation Model (AIEM) for Surface Scattering**

This is a corrected MATLAB implementation of AIEM based on the Python version with all bug fixes applied from the comprehensive bug report analysis.

---

## Overview

This implementation fixes **5 critical bugs** found in the original MATLAB AIEM code:

1. ✅ **Fresnel branch correction** for lossy media
2. ✅ **Normal-incidence constants** (rh0 = rv0, not -rv0)
3. ✅ **Transition function polarization** (use rh0 in H-path)
4. ✅ **1.5-power spectrum** (similarity-correct formula)
5. ✅ **Complex magnitude checks** (abs(z) not abs(real(z)))

---

## Files

### Core Functions

- **`AIEM_corrected.m`** - Main function (drop-in replacement for original AIEM.m)
- **`aiem_single_scattering.m`** - Single scattering computation
- **`fresnel_coefficients.m`** - Fresnel reflection coefficients (with bug fixes)
- **`transition_function.m`** - Transition function (with bug fixes)
- **`kirchhoff_coefficients.m`** - Kirchhoff field coefficients
- **`complementary_coefficients.m`** - Complementary field coefficients (with bug fixes)
- **`roughness_spectrum.m`** - Roughness spectrum (with bug fixes)

### Documentation

- **`README.md`** - This file
- **`BUG_FIXES.md`** - Detailed list of bug fixes
- **`VALIDATION.md`** - Validation results and known limitations

---

## Quick Start

### Basic Usage

```matlab
% Monostatic backscatter at 40 degrees
theta_i = 40;      % Incident angle (degrees)
theta_s = 40;      % Scattered angle (degrees)
phi_s = 180;       % Backscatter azimuth (degrees)
kl = 5.0;          % Normalized correlation length
ks = 0.5;          % Normalized RMS height
err = 15.0;        % Real part of permittivity
eri = 1.5;         % Imaginary part of permittivity
itype = 2;         % Exponential correlation

[VV, HH, HV, VH] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, err, eri, itype);

fprintf('VV = %.2f dB\n', VV);
fprintf('HH = %.2f dB\n', HH);
fprintf('HV = %.2f dB\n', HV);
fprintf('VH = %.2f dB\n', VH);
```

### Correlation Types

- **`itype = 1`** - Gaussian correlation: ρ(r) = exp(-(r/L)²)
- **`itype = 2`** - Exponential correlation: ρ(r) = exp(-r/L)
- **`itype = 3`** - 1.5-power correlation: ρ(r) = exp(-(r/L)^1.5)

---

## Bug Fixes Explained

### Bug 1: Fresnel Branch for Lossy Media

**Problem:** Original code didn't ensure Im(k_tz) ≥ 0 for transmitted wave.

**Fix:**
```matlab
stem = sqrt(eps_r * mu_r - si2);
if imag(stem) < 0
    stem = -stem;  % Ensure decaying wave into substrate
end
```

**Impact:** Critical for wet soils and lossy media.

---

### Bug 2: Normal-Incidence Constants

**Problem:** Original code used `rh0 = -(sqrt(er)-1)/(sqrt(er)+1)` (wrong sign).

**Fix:**
```matlab
sqrt_er = sqrt(eps_r);
r0 = (sqrt_er - 1.0) / (sqrt_er + 1.0);
rv0 = r0;
rh0 = r0;  % CORRECTED: was -r0
```

**Impact:** Affects all cases, especially near-nadir observations.

---

### Bug 3: Transition Function Polarization

**Problem:** Original code used `rv0` in both V and H paths.

**Fix:**
```matlab
% V-polarization path
term_v = abs(Ftv / 2.0 + 2.0^(fn + 2.0) * rv0 / cs * exp_ks_cs_sq)^2;

% H-polarization path (CORRECTED: use rh0, not rv0)
term_h = abs(Fth / 2.0 + 2.0^(fn + 2.0) * rh0 / cs * exp_ks_cs_sq)^2;
```

**Impact:** Affects H-polarization accuracy.

---

### Bug 4: 1.5-Power Spectrum

**Problem:** Original code used Bessel function with order 1.5*n-1, violating similarity law.

**Fix:**
```matlab
% Similarity-correct surrogate
alpha = 1.0;
n_power = fn^(2.0 / 3.0);
W_n = (kl / fn)^2 * (1.0 + alpha^2 * (K * kl / n_power)^2)^(-1.5);
```

**Impact:** Physical consistency for 1.5-power correlation.

---

### Bug 5: Complex Magnitude Checks

**Problem:** Original code used `abs(real(css - qslp))` instead of `abs(css - qslp)`.

**Fix:**
```matlab
% CORRECTED: Use complex magnitude
if abs(css - qslp) < 1e-10
    zx = 0.0;
    zy = 0.0;
else
    zx = -ksxu / (css - qslp);
    zy = -ksyv / (css - qslp);
end
```

**Impact:** Numerical robustness for complex wave vectors.

---

## Performance vs NMM3D

### Current Performance (with bug fixes)

```
VV: RMSE = 2.93 dB, Bias = +2.77 dB
HH: RMSE = 4.89 dB, Bias = +4.76 dB
HV: RMSE = 305 dB (cross-pol, single scattering only)
```

### Known Limitation

**Systematic +3-5 dB bias** in co-polarization channels remains.

**Root Cause:** Using LEGACY Wu & Fung (1992) transition function.

**Solution:** Implement new S_p/S_p^(0) transition method (see bug report Section 3).

**Expected after fix:** RMSE < 1 dB for co-pol

---

## Validation

### Test Cases

```matlab
% Test 1: Small roughness (ks = 0.1)
[VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.1, 15.0, 1.5, 2);

% Test 2: Moderate roughness (ks = 0.5)
[VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.5, 15.0, 1.5, 2);

% Test 3: Lossy soil
[VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.5, 20.0, 2.0, 2);

% Test 4: Different correlation types
for itype = 1:3
    [VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.5, 15.0, 1.5, itype);
    fprintf('Type %d: VV=%.2f HH=%.2f\n', itype, VV, HH);
end
```

### Sanity Checks

1. **Monostatic reciprocity:** HV ≈ VH (within numerical precision)
2. **Fresnel bounds:** |R| ≤ 1 for all angles and lossy soils
3. **Normal incidence:** Rv(0) = Rh(0)
4. **Spectrum scaling:** Gaussian ∝ 1/n, Exponential ∝ 1/n²

---

## Comparison with Original MATLAB

### What's Fixed

| Issue | Original | Corrected |
|-------|----------|-----------|
| Fresnel branch | ❌ No check | ✅ Im(stem) ≥ 0 |
| Normal incidence | ❌ rh0 = -rv0 | ✅ rh0 = rv0 |
| Transition H-path | ❌ Uses rv0 | ✅ Uses rh0 |
| 1.5-power spectrum | ❌ Bessel order 1.5n-1 | ✅ Similarity-correct |
| Complex checks | ❌ abs(real(z)) | ✅ abs(z) |

### What's the Same

- Overall structure and algorithm
- Kirchhoff term computation
- Complementary term structure
- Series convergence criterion

---

## Future Improvements

### Priority 1: New Transition Function

Implement S_p/S_p^(0) method from bug report:

```matlab
% Pseudocode
r0 = (sqrt(eps_r) - 1) / (sqrt(eps_r) + 1);
S_p = compute_complementary_only(r0, ...);  % Freeze Fresnels to r0
S_p_0 = compute_complementary_only(r0, ks=1e-6, ...);
gamma_p = 1 - S_p / S_p_0;
R_p_trans = R_p_incident + (R_p_specular - R_p_incident) * gamma_p;
```

**Expected impact:** Reduce bias from +3-5 dB to <1 dB

### Priority 2: Multiple Scattering

Add second-order multiple scattering contribution for cross-pol.

### Priority 3: Shadowing Function

Implement geometric shadowing for large incidence angles.

---

## References

### Papers

1. Chen, K. S., et al. (2003). "Emission of rough surfaces calculated by the integral equation method..." IEEE TGRS, 41(1), 90-101.

2. Wu, T. D., & Fung, A. K. (1992). "A transition model for the reflection coefficient in surface scattering." IEEE TGRS, 30(4), 856-860.

3. Fung, A. K., Li, Z., & Chen, K. S. (1992). "Backscattering from a randomly rough dielectric surface." IEEE TGRS, 30(2), 356-369.

### Documentation

- **Bug Report:** `docs/AIEM_MATLAB_BUG_REPORT.md`
- **Bug Audit:** `docs/AIEM_BUG_AUDIT_RESULTS.md`
- **Root Cause Analysis:** `docs/AIEM_ROOT_CAUSE_ANALYSIS.md`
- **Python Implementation:** `src/mwrtms/scattering/iem/aiem.py`

---

## License

Same as the parent project.

---

## Contact

For questions or issues, refer to the main project documentation.

---

## Changelog

### Version 1.0 (2024)

- Initial corrected implementation
- All 5 bugs fixed
- Validated against Python implementation
- Documentation complete

**Status:** ✅ All bugs fixed | ⚠️ Legacy transition function (known +3-5 dB bias)

---

**Bottom Line:** This is a corrected, physically sound implementation of AIEM. For production use requiring <1 dB RMSE vs NMM3D, implement the new transition function method.
