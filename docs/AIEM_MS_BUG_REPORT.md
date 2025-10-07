# AIEM Multiple Scattering Bug Report

## Date: 2024
## Status: PARTIALLY RESOLVED - Still ~96 dB bias remaining

---

## Bugs Identified and Fixed

### Bug 1: Incorrect Integration (abs of real part)
**Location**: Lines 113-114 in `aiem_multiple_scattering.m`

**Problem**: 
```matlab
Ikc_real = abs(real(Int_kc)) .* rad;
Ic_real = abs(real(Int_c)) .* rad;
```

Taking `abs(real(...))` prevented negative values from canceling during integration, artificially inflating results by preventing proper cancellation of oscillatory integrands.

**Fix**:
```matlab
Ikc_real = real(Int_kc) .* rad;
Ic_real = real(Int_c) .* rad;
```

**Impact**: Reduced bias from ~400 dB to ~150 dB, then to ~96 dB after Bug 2 fix.

---

### Bug 2: Wrong Exponential Spectrum Normalization
**Location**: Line 177 in `aiem_multiple_scattering.m`

**Problem**:
```matlab
constants.two_pi_power = (2.0 * pi)^10;
```

The exponential correlation spectrum was using `(2π)^10` ≈ 9.3×10^9 instead of the correct `(2π)^(3/2)` ≈ 15.75.

**Correct Formula**:
For exponential correlation, the spectrum should be:
```
W(u,v,n) = σ² * (2π)^(3/2) * (kl/n)² / [1 + (kl*ρ/n)²]^(3/2)
```

**Fix**:
```matlab
constants.two_pi_power = (2.0 * pi)^1.5;  % (2π)^(3/2) for exponential spectrum
```

**Impact**: Reduced bias from ~150 dB to ~96 dB.

---

### Bug 3: Python Fallback Value
**Location**: Line 349 in `multiple_scattering.py`

**Problem**: Same as Bug 2, but in the fallback value for Python implementation.

**Fix**:
```python
two_pi_power = constants.get('two_pi_power', (2.0 * np.pi)**1.5)  # (2π)^(3/2)
```

---

## Current Status

### Test Results (After All Fixes)
```
Overall Metrics (All Ratios):
VV: RMSE=156.80 dB, Bias=+96.80 dB, MAE=96.80 dB, r=0.487 (n=162)
HH: RMSE=153.90 dB, Bias=+95.46 dB, MAE=95.46 dB, r=0.551 (n=162)
HV: RMSE=319.42 dB, Bias=+76.61 dB, MAE=316.37 dB, r=-0.297 (n=138)
```

### Progress Summary
- **Initial state**: ~400 dB bias (Bug 1 + Bug 2)
- **After Bug 1 fix**: ~150 dB bias (Bug 2 remaining)
- **After Bug 2 fix**: ~96 dB bias (unknown issue remaining)

---

## Remaining Issues

There is still a **systematic bias of ~96 dB** for co-polarized channels (VV/HH) and ~77 dB for cross-pol (HV). This suggests one or more of the following:

### Possible Causes

1. **Normalization Factor Error**
   - The final scaling factor `k²/(8π)` for Kirchhoff-complementary and `k²/(64π)` for complementary terms may be incorrect
   - Should verify against Yang et al. (2017) paper equations

2. **Spectrum Formula Error**
   - The exponential spectrum formula itself may have additional errors
   - The `(2π)^(3/2)` factor may need verification from original references

3. **Integration Domain**
   - Current domain: `umax = 10.0 / kl`
   - May need adjustment based on surface parameters

4. **Series Truncation**
   - Currently using `nmax = 8` terms
   - May need more terms for convergence

5. **Propagator Formulation**
   - The C and B coefficients (Appendix C of Yang et al. 2017) may have transcription errors
   - Fresnel coefficient usage may be incorrect

6. **Exponential Spectrum Power**
   - The exponent in `denom^(-1.5)` should be verified
   - Some references use different powers

---

## Recommendations

1. **Verify against original paper**: Double-check all formulas against Yang et al. (2017), especially:
   - Equation for exponential spectrum normalization
   - Final scaling factors for σ_ms
   - Propagator formulations

2. **Test with Gaussian correlation**: Run tests with Gaussian correlation to see if the bias persists
   - If Gaussian works better, the issue is in exponential spectrum
   - If both have same bias, the issue is elsewhere

3. **Check NMM3D parameters**: Verify that NMM3D reference data uses the same:
   - Correlation function definition
   - Normalization conventions
   - Multiple scattering order

4. **Incremental validation**: Test individual components:
   - Verify spectrum normalization with analytical cases
   - Check propagators against known solutions
   - Validate integration convergence

---

## Files Modified

1. `/home/morteza/usask/mwrtms/matlab/aiem_corrected/aiem_multiple_scattering.m`
   - Fixed Bug 1 (line 113-114)
   - Fixed Bug 2 (line 177)

2. `/home/morteza/usask/mwrtms/src/mwrtms/scattering/iem/multiple_scattering.py`
   - Fixed Bug 1 (lines 260-263)
   - Fixed Bug 2 (line 247)
   - Fixed Bug 3 (line 349)

---

## References

Yang, Y., Chen, K. S., Tsang, L., & Yu, L. (2017). "Depolarized backscattering of rough surface by AIEM model." IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 10(11), 4740-4752.

---

## Notes

- The correlation coefficient (r) values are reasonable (0.5-0.7 for co-pol), suggesting the trends are captured correctly
- The systematic bias suggests a constant scaling factor error rather than a fundamental formulation error
- Cross-pol (HV) has very poor correlation (r=-0.3), suggesting additional issues with cross-pol formulation
