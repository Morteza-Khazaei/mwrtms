# HH Band Bug Fix: I2EM Transition in AIEM

## The Problem

When `use_i2em_transition=True` was enabled:
- ✅ **VV improved**: RMSE dropped from 5.7 dB to 2.1 dB
- ❌ **HH got WORSE**: RMSE increased from 4.0 dB to 6.2 dB  
- ❌ **HV completely broken**: 306 dB error

## Root Cause

**Line 193 in aiem.py** was using I2EM's negative `rh0` convention:

```python
if self._use_i2em_transition:
    Tfv, Tfh = compute_i2em_transition_function(...)
    rh0_trans = -rv0  # ← BUG: This is wrong for AIEM!
```

### Why This Caused the Problem

1. **I2EM uses `rh0 = -rv0`** (negative nadir reflection for HH)
2. **AIEM uses `rh0 = +rv0`** (positive nadir reflection for HH)
3. **AIEM's complementary terms expect positive `rh0`**

When we used negative `rh0` in AIEM:
- `Rhtran` stayed very negative (-0.643 instead of +0.505)
- This increased `|fhh|²` by 1.62x (2.09 dB)
- The complementary terms still used `Rhi`, creating a mismatch
- **Result**: HH backscatter increased instead of decreased!

### Numerical Example (40° incidence, ε=15+3j)

```
Fresnel coefficients:
  Rhi = -0.670 (negative at incident angle)
  rv0 = +0.594 (positive at nadir)

AIEM original:
  rh0 = +0.594
  Rhtran = Rhi + (rh0 - Rhi) * Tfh
         = -0.670 + (+0.594 - (-0.670)) * 0.93
         = +0.505  ← Becomes positive!

I2EM transition (WRONG):
  rh0 = -0.594  ← Negative!
  Rhtran = Rhi + (rh0 - Rhi) * Tf
         = -0.670 + (-0.594 - (-0.670)) * 0.35
         = -0.643  ← Stays very negative!

Impact:
  |Rhtran|² ratio = 0.414 / 0.256 = 1.62 (2.09 dB increase!)
```

## The Fix

**Use I2EM's transition factor but keep AIEM's positive `rh0` convention:**

```python
if self._use_i2em_transition:
    Tfv, Tfh = compute_i2em_transition_function(...)
    # Use I2EM transition factor but keep AIEM's rh0 convention
    # (AIEM's complementary terms expect positive rh0)
    rh0_trans = rh0  # ← FIXED: Use positive rh0
```

This gives us:
- ✅ I2EM's better transition factor (Tf = 0.35 instead of 0.81-0.93)
- ✅ AIEM's rh0 convention (positive, compatible with complementary terms)
- ✅ Best of both worlds!

## Expected Results After Fix

### VV Band (Should Stay Good)
- Before fix: RMSE = 2.1 dB ✓
- After fix: RMSE ≈ 2.1 dB (no change expected)

### HH Band (Should Improve!)
- Before fix: RMSE = 6.2 dB ❌
- After fix: RMSE ≈ 2-3 dB ✓ (similar to original AIEM or better)

### HV Band (Should Fix)
- Before fix: RMSE = 306 dB ❌ (completely broken)
- After fix: RMSE ≈ 5-10 dB ✓ (back to normal)

## Why This Works

The key insight is that **I2EM and AIEM use different conventions**:

| Aspect | I2EM | AIEM |
|--------|------|------|
| **Transition factor** | Tf = 0.35 (better!) | Tfv = 0.81, Tfh = 0.93 (too large) |
| **rh0 convention** | -rv0 (negative) | +rv0 (positive) |
| **Complementary terms** | Simple | Complex (8 branches) |

We can **mix and match**:
- Take I2EM's better transition factor ✓
- Keep AIEM's rh0 convention ✓
- Use AIEM's complementary terms ✓

This hybrid approach gives optimal performance!

## Validation

To confirm the fix works, re-run the notebook tests and check:

1. **VV should stay good**: RMSE ≈ 1.5-2.0 dB
2. **HH should improve**: RMSE ≈ 2-3 dB (not 6 dB!)
3. **HV should work**: RMSE ≈ 5-10 dB (not 306 dB!)

## Lesson Learned

**Don't blindly copy conventions between models!**

I2EM and AIEM have different internal conventions. When borrowing the transition function from I2EM, we need to:
1. ✅ Take the transition **algorithm** (how Tf is computed)
2. ❌ Don't take the **sign conventions** (rh0 = -rv0 vs +rv0)
3. ✅ Keep compatibility with the rest of AIEM's code

The fix is simple but critical: use I2EM's Tf but AIEM's rh0!
