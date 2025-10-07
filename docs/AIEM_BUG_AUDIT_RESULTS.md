# AIEM Bug Audit Results

This document compares the MATLAB AIEM bug report against our Python implementation to identify which bugs have been addressed and which require fixes.

## Summary

| Bug # | Issue | Status | Action Required |
|-------|-------|--------|-----------------|
| 1 | Specular half-angle formula | ✅ CORRECT | None - already correct |
| 2 | Fresnel branch for lossy media | ❌ MISSING | **FIX REQUIRED** |
| 3 | Normal-incidence constants | ❌ WRONG | **FIX REQUIRED** |
| 4 | Transition function | ❌ LEGACY | **FIX REQUIRED** |
| 5 | 1.5-power spectrum | ❌ WRONG | **FIX REQUIRED** |
| 6 | Complex near-singularity guards | ⚠️ PARTIAL | Minor improvement needed |
| 7 | Bessel symmetry | ✅ CORRECT | None - using scipy |

---

## Detailed Analysis

### Bug 1: Specular Half-Angle Formula ✅ CORRECT

**Report says:** The formula `csl = sqrt((1 + cs*css - si*sis*csfs)/2)` is correct.

**Our code** (`fresnel_utils.py` line 77):
```python
csl = np.sqrt(1.0 + cs * css - si * sis * csfs) / np.sqrt(2.0)
```

**Status:** ✅ **CORRECT** - No action needed.

---

### Bug 2: Fresnel Branch for Lossy Media ❌ MISSING

**Report says:** 
> Use k_tz/k = sqrt(ε_r - sin²θ) with the **decaying branch**: if Im{k_tz} < 0, flip the sign.

**Our code** (`fresnel_utils.py` lines 54-56):
```python
# Transmitted wave vector component
stem = np.sqrt(eps_r * mu_r - si2)
```

**Problem:** No branch correction for lossy media. When `Im(stem) < 0`, we need to flip the sign to ensure the transmitted wave decays into the lower half-space.

**Status:** ❌ **FIX REQUIRED**

**Fix:**
```python
# Transmitted wave vector component with proper branch
stem = np.sqrt(eps_r * mu_r - si2)
# Ensure decaying branch: Im(stem) >= 0
if np.imag(stem) < 0:
    stem = -stem
```

---

### Bug 3: Normal-Incidence Constants ❌ WRONG

**Report says:**
> Set r_0 = (1 - sqrt(ε_r))/(1 + sqrt(ε_r)) and use **the same r_0** for both pols.
> **Bug in MATLAB:** `rv0 = (sqrt(er)-1)/(sqrt(er)+1)`, `rh0 = -(sqrt(er)-1)/(sqrt(er)+1)`

**Our code** (`fresnel_utils.py` lines 135-138):
```python
sqrt_er = np.sqrt(eps_r)
r0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
rv0 = r0
rh0 = r0  # Same as rv0 at normal incidence (CORRECTED from MATLAB bug)
```

**Good!** But check `transition.py` (lines 60-61):
```python
rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0)
rh0 = -rv0  # ❌ WRONG!
```

**Problem:** The transition function still uses the MATLAB bug with opposite signs!

**Status:** ❌ **FIX REQUIRED in transition.py**

**Fix:** Change line 61 in `transition.py` from:
```python
rh0 = -rv0
```
to:
```python
rh0 = rv0  # Same as rv0 at normal incidence
```

---

### Bug 4: Transition Function ❌ LEGACY

**Report says:**
> Implement R_p^(T) = R_p(θ_i) + [R_p(θ_sp) - R_p(θ_i)](1 - S_p/S_p^(0))
> 
> If legacy block is kept:
> - Fix the pol typo: in H-path use `rh0` not `rv0`
> - Make Gaussian factors consistent: multiply by exp(-(ks·cs)²) in all places

**Our code** (`transition.py`):
Uses the legacy Wu & Fung formulation with the shadowing terms St/St0.

**Problems identified:**
1. Line 61: Uses `rh0 = -rv0` (wrong sign, see Bug 3)
2. Lines 88-89: Uses `rv0` in both V and H paths:
   ```python
   term_v = np.abs(Ftv / 2.0 + (2.0 ** (fn + 2.0)) * rv0 / cs * exp_ks_cs_sq) ** 2
   term_h = np.abs(Fth / 2.0 + (2.0 ** (fn + 2.0)) * rv0 / cs * exp_ks_cs_sq) ** 2
   ```
   Should use `rh0` in the H-path.

3. Gaussian factor: Uses `exp_ks_cs_sq = np.exp(-ks_cs_sq)` which is correct, but need to verify consistency.

**Status:** ❌ **FIX REQUIRED**

**Recommendation:** The report suggests implementing the new S_p/S_p^(0) method, but if keeping legacy:
- Fix `rh0 = rv0` (not `-rv0`)
- Use `rh0` in H-path (line 89)

---

### Bug 5: 1.5-Power Spectrum ❌ WRONG

**Report says:**
> **Bug in MATLAB:** current code uses modified-Bessel form with **order 1.5*n-1**. That violates the similarity law.
>
> **Fix:** Remove the Bessel-K form with order 1.5*n-1. Replace with either:
> - Numerical Hankel transform (exact)
> - Similarity-correct surrogate: W^(n)(K) ≈ (L/n)² (1 + α²(KL/n^(2/3))²)^(-1.5)

**Our code** (`spectrum_aiem.py` lines 109-123):
```python
elif corr_type in ('3', 'powerlaw', 'power', '1.5-power', 'xpower'):
    # 1.5-power (transformed exponential) correlation function
    if np.isclose(K, 0.0):
        spectrum = kl2 / (3.0 * fn - 2.0)
    else:
        e = power_exponent * fn - 1.0
        y = power_exponent * fn
        m = power_exponent * fn - 1.0  # ❌ Order depends on n!
        
        # Use log-domain to avoid overflow
        try:
            log_gamma = np.log(gamma(y))
            log_bessel = np.log(kv(m, K))  # ❌ Bessel order = 1.5*n - 1
            ...
```

**Problem:** The Bessel order `m = 1.5*n - 1` depends on `n`, which violates the similarity law. The spectrum should have the form:
```
W^(n)(K) = L² * n^(-4/3) * Φ(K*L*n^(-2/3))
```
where Φ is independent of n.

**Status:** ❌ **FIX REQUIRED**

**Fix:** Implement the similarity-correct surrogate:
```python
# Similarity-correct surrogate for 1.5-power spectrum
alpha = 1.0  # Tuning parameter to match curvature at K=0
spectrum = (kl / fn) ** 2 * (1.0 + alpha**2 * (K * kl / fn**(2.0/3.0))**2) ** (-1.5)
```

Or implement numerical Hankel transform (more complex but exact).

---

### Bug 6: Complex Near-Singularity Guards ⚠️ PARTIAL

**Report says:**
> Replace tests like `abs(real(css-qslp))<tol` by `abs(css-qslp)<tol` (complex magnitude).

**Our code:** Need to check complementary field coefficient functions.

**Status:** ⚠️ Need to audit complementary.py for this issue.

---

### Bug 7: Bessel Symmetry ✅ CORRECT

**Report says:**
> Prefer K_{|ν|}(x) numerically (`besselk(abs(nu),x)`), since K_{-ν} = K_ν.

**Our code:** Uses `scipy.special.kv` which handles this correctly.

**Status:** ✅ **CORRECT** - No action needed.

---

## Priority Fixes

### HIGH PRIORITY (Affects accuracy for lossy soils and all cases)

1. **Bug 2: Fresnel branch correction** - Critical for lossy media
2. **Bug 3: Normal-incidence constants in transition.py** - Affects all cases
3. **Bug 4: Transition function pol typo** - Affects H-pol accuracy

### MEDIUM PRIORITY (Affects specific correlation functions)

4. **Bug 5: 1.5-power spectrum** - Only affects users of this correlation type

### LOW PRIORITY (Minor improvements)

5. **Bug 6: Complex magnitude checks** - Numerical robustness

---

## Implementation Plan

1. Fix `fresnel_utils.py`:
   - Add branch correction for `stem` in both `compute_fresnel_incident` and `compute_fresnel_specular`

2. Fix `transition.py`:
   - Change `rh0 = -rv0` to `rh0 = rv0`
   - Change line 89 to use `rh0` instead of `rv0`

3. Fix `spectrum_aiem.py`:
   - Replace 1.5-power spectrum with similarity-correct formula

4. Audit `complementary.py`:
   - Check for `abs(real(...))` patterns and replace with `abs(...)`

5. Add tests to verify:
   - Monostatic HV = VH
   - Lossy soil behavior
   - Small roughness limits
   - Spectrum scaling laws

---

## Testing Strategy

After fixes, validate:

1. **Fresnel coefficients:** |R_p(θ)| ≤ 1 for all angles and lossy soils
2. **Normal incidence:** R_h(0) = R_v(0) = r_0
3. **Monostatic reciprocity:** σ_hv = σ_vh
4. **Small roughness:** γ_p finite, S_p^(0) stable at ks → 0
5. **Spectrum scaling:**
   - Gaussian: width ∝ √n
   - Exponential: amplitude ∝ 1/n²
   - 1.5-power: n^(4/3) * W^(n) vs K*L*n^(-2/3) should collapse

---

## Conclusion

Our Python implementation has **4 critical bugs** that need immediate fixing:

1. Missing Fresnel branch correction for lossy media
2. Wrong sign for rh0 in transition function
3. Wrong polarization constant in H-path of transition
4. Wrong 1.5-power spectrum formula

The specular angle formula and overall structure are correct. After these fixes, the implementation will be physically sound and ready for production use.
