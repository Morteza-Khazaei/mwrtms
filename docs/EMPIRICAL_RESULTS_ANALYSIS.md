# Empirical Results: AIEM vs I2EM vs NMM3D

## Actual Test Results (162 Cases)

### VV Polarization
```
AIEM:  RMSE = 14.81 dB  |  Bias = -14.54 dB  |  Corr = 0.942
I2EM:  RMSE = 16.70 dB  |  Bias = -16.47 dB  |  Corr = 0.963
```
**Winner: AIEM** (1.89 dB lower RMSE)

### HH Polarization
```
AIEM:  RMSE = 16.09 dB  |  Bias = -15.96 dB  |  Corr = 0.988
I2EM:  RMSE = 19.32 dB  |  Bias = -19.21 dB  |  Corr = 0.983
```
**Winner: AIEM** (3.23 dB lower RMSE)

---

## Critical Findings

### 1. AIEM Actually Performs BETTER Than I2EM!

This contradicts the documented behavior. AIEM has:
- **Lower RMSE** for both VV and HH
- **Less bias** (closer to zero)
- **Good correlation** (0.94-0.99)

### 2. Both Models Show NEGATIVE Bias (Under-Prediction)

```
Expected:  AIEM has +3-5 dB positive bias (over-prediction)
Actual:    AIEM has -14 to -16 dB negative bias (under-prediction)
           I2EM has -16 to -19 dB negative bias (under-prediction)
```

This is the **opposite** of what was documented!

### 3. The Bias is HUGE

Both models under-predict by 14-19 dB. This suggests:
- **Systematic implementation error** in both models
- **Wrong normalization factor**
- **Missing physics** (e.g., multiple scattering not included)
- **Data format mismatch** between models and NMM3D

---

## What This Means

### The Original Analysis Was Wrong

The documentation claimed:
- ❌ "AIEM has +3-5 dB positive bias"
- ❌ "I2EM performs better for co-pol"
- ❌ "AIEM's transition function causes over-prediction"

**Reality**:
- ✅ AIEM has -14 to -16 dB negative bias
- ✅ AIEM performs BETTER than I2EM
- ✅ Both models severely UNDER-predict

### The Real Problem

The issue is NOT that AIEM over-predicts. The issue is that **BOTH models under-predict by ~15-19 dB**!

Possible causes:
1. **Normalization error**: Missing factor of ~30-100 (10-20 dB)
2. **Multiple scattering**: Not included (could add 5-10 dB)
3. **Shadowing**: Not properly applied
4. **Units mismatch**: dB vs linear, or area normalization

---

## Why the Confusion?

The documentation may have been based on:
1. **Different test conditions** (different roughness range)
2. **Older AIEM version** with different bugs
3. **Misinterpretation** of earlier results
4. **Theoretical expectations** rather than empirical tests

---

## Next Steps

### 1. Investigate the 15-19 dB Under-Prediction

This is the REAL problem, not the transition function!

Check:
- Normalization factors in both models
- Multiple scattering contribution
- Shadowing function application
- Unit conversions

### 2. Test with Multiple Scattering Enabled

```python
model_aiem = AIEMModel(
    ...,
    include_multiple_scattering=True  # Add MS contribution
)
```

This might reduce the bias significantly.

### 3. Compare Against Known-Good Implementation

Test against MATLAB reference code to see if the issue is in:
- The Python implementation
- The model formulation itself
- The test setup

### 4. Check Specific Roughness Ranges

The documented +3-5 dB bias might occur only for:
- Specific ks ranges
- Specific kL ranges
- Specific permittivity values

---

## Conclusion

**The empirical test reveals**:

1. ✅ **AIEM performs BETTER than I2EM** (not worse!)
   - VV: 1.89 dB lower RMSE
   - HH: 3.23 dB lower RMSE

2. ❌ **Both models severely under-predict** (~15-19 dB bias)
   - This is the REAL problem to solve
   - Not the transition function issue

3. ⚠️ **The documented behavior is incorrect**
   - No +3-5 dB positive bias found
   - AIEM doesn't over-predict
   - I2EM doesn't perform better

**The I2EM transition implementation was based on incorrect assumptions!**

The real issue to investigate is: **Why do both models under-predict by 15-19 dB?**

This could be:
- Missing multiple scattering (5-10 dB)
- Wrong normalization (10-20 dB)
- Shadowing not applied correctly
- Some other systematic error

---

## Recommendation

**STOP** focusing on the transition function. Instead:

1. **Find the 15-19 dB missing power**
   - Check normalization factors
   - Enable multiple scattering
   - Verify shadowing functions

2. **Compare with MATLAB reference**
   - Run same test cases in MATLAB AIEM
   - See if MATLAB also under-predicts

3. **Test different roughness ranges**
   - The +3-5 dB bias might exist for specific conditions
   - Current test uses full NMM3D dataset

4. **Verify NMM3D data format**
   - Ensure we're reading it correctly
   - Check units and conventions

The transition function work was premature optimization based on incorrect problem diagnosis!
