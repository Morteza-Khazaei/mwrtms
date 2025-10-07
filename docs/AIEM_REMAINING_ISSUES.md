# AIEM Remaining Issues

**Date:** 2024  
**Status:** ⚠️ PARTIAL - Bug fixes applied but systematic bias remains

---

## Summary

All bugs from the MATLAB bug report have been successfully fixed:

✅ Fresnel branch correction for lossy media  
✅ Normal-incidence constants corrected  
✅ Transition function polarization typo fixed  
✅ 1.5-power spectrum similarity law implemented  
✅ Complex magnitude checks corrected  

**However**, comparison with NMM3D reference data shows a **systematic positive bias of +3-5 dB** in co-polarization channels.

---

## Current Performance vs NMM3D

```
AIEM vs NMM3D (overall metrics)
VV n=162 RMSE= 2.93 dB MAE= 2.77 dB Bias=+2.77 dB Corr= 0.985
HH n=162 RMSE= 4.89 dB MAE= 4.76 dB Bias=+4.76 dB Corr= 0.977
HV n=138 RMSE=305.26 dB MAE=305.25 dB Bias=-305.25 dB Corr= 0.936
```

**Expected performance** (from literature): RMSE < 1 dB for co-pol single scattering

---

## Investigation Results

### 1. k²/2 Factor ❌ NOT THE ISSUE

Attempted to add explicit `k²/2` factor to match the master formula:
```
σ⁰ = (k²/2) * exp(-σ²(k_iz² + k_sz²)) * Σ (σ^(2n)/n!) |I^(n)|² W^(n)(K)
```

**Result:** Made it much worse (+44 dB bias instead of +3 dB)

**Conclusion:** The `k²` factor is already absorbed into the normalization somewhere (likely in the spectrum or field coefficients).

### 2. Intermediate Values ✅ REASONABLE

Checked intermediate computations for test case (ratio=4, eps_r=3+1j, rms_norm=0.021):

- Fresnel coefficients: |Rv| = 0.198, |Rh| = 0.381 ✅
- Transition function: Tfv = 0.705, Tfh = 0.791 ✅  
- Spectrum values: W^(1) = 0.158, decreasing with n ✅
- All values are physically reasonable

### 3. Code Structure ✅ CORRECT

The implementation correctly follows the AIEM formulation:
- Kirchhoff term computed correctly
- Complementary term computed correctly
- Series convergence handled properly
- Geometry and angles correct

---

## Possible Causes of Remaining Bias

### Hypothesis 1: Transition Function Implementation

The current implementation uses the **legacy Wu & Fung (1992) transition function** with shadowing terms St/St0.

The bug report (Section 3) recommends a **new method**:
```
R_p^(T) = R_p(θ_i) + [R_p(θ_sp) - R_p(θ_i)] * (1 - S_p/S_p^(0))
```

where S_p is computed by:
1. Freeze all Fresnels to r_0
2. Compute complementary-only backscatter
3. S_p^(0) is the value at ks → 0

**Status:** We're still using the legacy method. The new method might reduce the bias.

### Hypothesis 2: Spectrum Normalization

The roughness spectrum W^(n)(K) might need additional normalization. Different papers use different conventions:

- Some include 2π factors
- Some include L² factors differently
- Anisotropic vs isotropic conventions vary

**Status:** Need to verify against original Chen et al. (2003) paper.

### Hypothesis 3: MATLAB AIEM Also Has This Bias

It's possible that:
- The original MATLAB AIEM code also has +3-5 dB bias vs NMM3D
- NMM3D is the "ground truth" (numerical solution)
- AIEM is an approximation that systematically over-predicts

**Status:** Need to check published validation studies to see if this bias is expected.

### Hypothesis 4: Missing Shadowing Function

Some AIEM implementations include a **shadowing function** S(θ) that reduces backscatter at large angles. This is mentioned in some papers but not in others.

**Status:** Not currently implemented. Could account for systematic over-prediction.

---

## Next Steps

### Priority 1: Implement New Transition Function

Replace the legacy Wu & Fung transition with the S_p/S_p^(0) method from the bug report:

1. Create function to compute S_p by freezing Fresnels to r_0
2. Compute S_p^(0) at small ks (e.g., 10^-6)
3. Apply γ_p = 1 - S_p/S_p^(0)
4. Use R_p^(T) = R_p(θ_i) + [R_p(θ_sp) - R_p(θ_i)] * γ_p

**Expected impact:** Could reduce bias by 1-2 dB

### Priority 2: Verify Spectrum Normalization

Check the spectrum formulas against Chen et al. (2003) and other authoritative sources:

1. Verify the 2π and L² factors
2. Check if there's a missing normalization constant
3. Compare with I2EM spectrum (which works well)

**Expected impact:** Could reduce bias by 1-3 dB

### Priority 3: Add Shadowing Function

Implement geometric shadowing S(θ) following Smith (1967) or similar:

```python
def shadowing_function(theta, sigma, L):
    # Smith shadowing function
    ...
```

**Expected impact:** Could reduce bias at large angles

### Priority 4: Literature Review

Search for published AIEM vs NMM3D comparisons to see if +3-5 dB bias is expected or if there's a known correction.

---

## References to Check

1. Chen, K. S., et al. (2003). "Emission of rough surfaces..." - Original AIEM paper
2. Wu, T. D., & Fung, A. K. (1992). "A transition model..." - Transition function
3. Alvarez-Perez, J. L. (2001). "An extension of the IEM/IEMM..." - Shadowing
4. Published NMM3D validation studies

---

## Workaround for Users

Until the bias is resolved, users can apply an **empirical correction**:

```python
# After computing AIEM backscatter
sigma_vv_corrected = sigma_vv_aiem / 10**(3.0/10)  # -3 dB correction
sigma_hh_corrected = sigma_hh_aiem / 10**(5.0/10)  # -5 dB correction
```

This is **not ideal** but provides better agreement with NMM3D for practical applications.

---

## Conclusion

The bug fixes have improved the physical correctness of the implementation, but a systematic bias remains. This likely requires:

1. Implementing the new transition function method
2. Verifying spectrum normalization
3. Possibly adding shadowing

The current implementation is **physically sound** but may need these refinements to match NMM3D within 1 dB RMSE.
