# Empirical Test Results: Reality Check

## Executive Summary

**The empirical test contradicts all previous assumptions!**

- ✅ **AIEM performs BETTER than I2EM** (not worse)
- ❌ **Both models severely UNDER-predict** (not over-predict)
- ⚠️ **The real problem is 15-19 dB missing power** (not transition function)

---

## Actual Test Results (162 NMM3D Cases)

### VV Polarization
```
AIEM:  RMSE = 14.81 dB  |  Bias = -14.54 dB  |  Corr = 0.942
I2EM:  RMSE = 16.70 dB  |  Bias = -16.47 dB  |  Corr = 0.963

Winner: AIEM (1.89 dB lower RMSE)
```

### HH Polarization
```
AIEM:  RMSE = 16.09 dB  |  Bias = -15.96 dB  |  Corr = 0.988
I2EM:  RMSE = 19.32 dB  |  Bias = -19.21 dB  |  Corr = 0.983

Winner: AIEM (3.23 dB lower RMSE)
```

---

## What We Expected vs What We Found

| Aspect | Expected (Documented) | Actual (Empirical) |
|--------|----------------------|-------------------|
| **AIEM Bias** | +3 to +5 dB (over-prediction) | **-14 to -16 dB (under-prediction)** |
| **I2EM Bias** | ~0 dB (accurate) | **-16 to -19 dB (under-prediction)** |
| **Better Model** | I2EM | **AIEM** |
| **RMSE Difference** | I2EM 2-3 dB better | **AIEM 2-3 dB better** |
| **Problem** | Transition too aggressive | **Both models missing 15-19 dB!** |

---

## Critical Findings

### 1. AIEM Performs Better Than I2EM

**Contrary to documentation**, AIEM has:
- Lower RMSE (14.81 vs 16.70 dB for VV)
- Less bias (-14.54 vs -16.47 dB for VV)
- Good correlation (0.94-0.99)

### 2. Both Models Severely Under-Predict

**Both models are missing 15-19 dB of backscatter power!**

This is NOT a transition function issue - it's a fundamental implementation problem affecting both models equally.

### 3. The Sign is Wrong

- Expected: Positive bias (over-prediction)
- Actual: **Negative bias (under-prediction)**

This suggests the problem is completely different from what was documented.

---

## Why the Confusion?

The previous analysis was based on:
1. **Theoretical expectations** rather than empirical tests
2. **Incomplete testing** or different test conditions
3. **Misinterpretation** of earlier results
4. **Documentation from different AIEM version**

**The I2EM transition implementation was solving the WRONG problem!**

---

## The Real Problem

### Missing ~15-19 dB of Power

Possible causes:

1. **Normalization Factor** (10-20 dB)
   - Wrong k² factor
   - Missing 2π or 4π
   - Area normalization issue

2. **Multiple Scattering** (5-10 dB)
   - Not included in test
   - Could add significant power

3. **Shadowing Function** (2-5 dB)
   - Not applied correctly
   - Wrong implementation

4. **Units/Conventions** (variable)
   - dB vs linear mismatch
   - Steradian vs square meter
   - Different reference systems

---

## What to Do Next

### 1. Enable Multiple Scattering

```python
model_aiem = AIEMModel(
    wave, geometry, surface,
    include_multiple_scattering=True  # This might add 5-10 dB
)
```

Test if this reduces the bias.

### 2. Check Normalization Factors

Review the final scaling in both models:
- AIEM: `0.5 * exp(...) * sum_val`
- I2EM: `k**2 / 2.0 * exp(...)`

Compare with reference implementations.

### 3. Verify Against MATLAB

Run the same test cases in MATLAB AIEM to see if:
- MATLAB also under-predicts
- The issue is in Python implementation
- The issue is in the model itself

### 4. Test Specific Roughness Ranges

The documented +3-5 dB bias might exist only for:
- Specific ks ranges (e.g., ks < 1.0)
- Specific kL ranges (e.g., kL > 10)
- Specific conditions not in full NMM3D dataset

---

## Implications

### The Transition Function Work Was Premature

All the analysis about:
- 2^(n+2) vs 2^(n+1)
- Transition factor being too large
- I2EM's better empirical tuning

**Was based on incorrect problem diagnosis!**

The real issue is that **both models are missing 15-19 dB**, not that one has a 3-5 dB bias relative to the other.

### AIEM is Actually Better

The empirical test shows AIEM:
- Has 2-3 dB lower RMSE than I2EM
- Is closer to NMM3D (less bias)
- Has excellent correlation

**AIEM doesn't need I2EM's transition - it's already better!**

### Both Models Need Fixing

The priority should be:
1. Find the missing 15-19 dB
2. Fix the fundamental implementation issue
3. THEN worry about fine-tuning (if needed)

---

## Conclusion

**The empirical test reveals the truth:**

1. ✅ **AIEM > I2EM** (contrary to documentation)
2. ❌ **Both under-predict by 15-19 dB** (the real problem)
3. ⚠️ **Transition function work was premature** (wrong problem)

**Next steps:**
- Enable multiple scattering and retest
- Check normalization factors
- Compare with MATLAB reference
- Find the missing 15-19 dB!

**The I2EM transition implementation, while technically correct, was solving a problem that doesn't exist in the empirical data!**

---

## Test Details

- **Dataset**: NMM3D_LUT_NRCS_40degree.dat (162 cases)
- **Frequency**: 5.4 GHz (C-band)
- **Incidence angle**: 40°
- **Correlation**: Exponential
- **Spectral terms**: 10 (fixed)
- **Multiple scattering**: Disabled
- **Success rate**: 100% (162/162 cases)

Both models succeeded on all test cases, allowing direct comparison.
