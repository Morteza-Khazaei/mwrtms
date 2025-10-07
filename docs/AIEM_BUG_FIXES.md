# AIEM Bug Fixes Based on MATLAB Bug Report

## Summary of Bugs Found and Fixed

### ✅ Bug #1: Specular Half-Angle Sign Error (CRITICAL)

**Location:** `src/mwrtms/scattering/iem/fresnel_utils.py` line 95

**Bug:** 
```python
csl = np.sqrt(1.0 + cs * css - si * sis * csfs) / np.sqrt(2.0)
```

**Should be:**
```python
csl = np.sqrt(1.0 - cs * css + si * sis * csfs) / np.sqrt(2.0)
```

**Explanation:** The formula for the specular half-angle uses cos(ψ/2) = sqrt((1 + cos(ψ))/2), where cos(ψ) = -cos(θ_i)cos(θ_s) + sin(θ_i)sin(θ_s)cos(φ_s). The MATLAB code had the wrong sign pattern.

---

### ✅ Bug #2: Fresnel Branch for Complex Permittivity

**Status:** ✅ **ALREADY CORRECT** in our implementation

**Location:** `src/mwrtms/scattering/iem/fresnel_utils.py`

Our code uses `np.sqrt()` which automatically handles the complex branch correctly for NumPy complex numbers. The principal square root in NumPy ensures Im(sqrt(z)) ≥ 0 for complex z, which is the physically correct branch for decay into the substrate.

---

### ✅ Bug #3: Normal-Incidence Constants (CRITICAL)

**Location:** `src/mwrtms/scattering/iem/fresnel_utils.py` lines 130-135

**Bug:**
```python
rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
rh0 = -rv0  # WRONG!
```

**Should be:**
```python
r0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
rv0 = r0
rh0 = r0  # Same for both polarizations at normal incidence
```

**Explanation:** At normal incidence, there is no distinction between H and V polarizations. Both should have the same value.

---

### ⚠️ Bug #4: Transition Function

**Status:** ⚠️ **NEEDS INVESTIGATION**

**Location:** `src/mwrtms/scattering/iem/transition.py`

The bug report recommends replacing the legacy transition function with:
```
R_p^(T) = R_p(θ_i) + (R_p(θ_sp) - R_p(θ_i)) * γ_p
where γ_p = 1 - S_p / S_p^(0)
```

Our current implementation follows the MATLAB approach. This needs careful review and testing before changing, as it affects all results.

---

### ✅ Bug #5: 1.5-Power Spectrum

**Status:** ✅ **NOT USED** in our implementation

We only implement Gaussian and Exponential correlation functions. The 1.5-power spectrum is not implemented, so this bug doesn't affect us.

---

### ✅ Bug #6: Complex Near-Singularity Guards

**Status:** ✅ **ALREADY CORRECT**

**Location:** Throughout the code

Our Python implementation uses `np.abs()` for complex magnitude checks, which is correct. Example from `geometry_utils.py`:
```python
if np.abs(css - qslp) < tolerance:
    # Handle singularity
```

---

### ✅ Bug #7: Bessel Symmetry

**Status:** ✅ **NOT APPLICABLE**

We don't use Bessel functions in our Gaussian/Exponential implementations.

---

## Priority Fixes Required

### HIGH PRIORITY (Affects Results):
1. ✅ **Fix specular half-angle** (Bug #1)
2. ✅ **Fix nadir coefficients** (Bug #3)

### MEDIUM PRIORITY (Review Needed):
3. ⚠️ **Review transition function** (Bug #4) - requires careful testing

### LOW PRIORITY (Already Correct):
4. ✅ Fresnel branch handling (Bug #2)
5. ✅ Singularity guards (Bug #6)
6. ✅ 1.5-power spectrum (Bug #5) - not implemented
7. ✅ Bessel symmetry (Bug #7) - not applicable

