# HH Failure at High Roughness with I2EM Transition

## The Problem

With `use_i2em_transition=True`, HH performance degrades catastrophically at higher roughness:

| kL Band | ks Range | HH RMSE | HH Correlation | Status |
|---------|----------|---------|----------------|--------|
| 4.0 | Low | 5.81 dB | 0.970 | Acceptable |
| 7.0 | Medium | 6.15 dB | **0.328** | Poor |
| 10.0 | High | 7.20 dB | **0.005** | Failed |
| 15.0 | Very High | 8.50 dB | **-0.212** | Completely Failed |

**The correlation drops from 0.97 to -0.21 as roughness increases!**

## Root Cause Analysis

### The I2EM Transition Factor is Too Small for HH

From the diagnostic, we know:
- I2EM uses **single Tf = 0.35** for both V and H
- AIEM uses **Tfv = 0.81, Tfh = 0.93** (different for each)

When we use I2EM's Tf = 0.35 for HH:
```
Rhtran = Rhi + (rh0 - Rhi) * 0.35
```

At **low roughness** (ks < 0.5):
- Rhi ≈ -0.67
- rh0 ≈ +0.59
- Rhtran ≈ -0.67 + 1.26 * 0.35 = -0.23 (acceptable)

At **high roughness** (ks > 1.0):
- Rhi becomes more negative
- The small Tf = 0.35 doesn't transition enough
- Rhtran stays too negative
- **The model diverges from NMM3D**

### Why AIEM's Original Tfh = 0.93 Works Better for HH

AIEM uses a **larger transition factor for HH** (0.93 vs 0.81 for VV):
```
Rhtran = Rhi + (rh0 - Rhi) * 0.93
       = -0.67 + 1.26 * 0.93
       = +0.50  ← Much more positive!
```

This stronger transition is **necessary for HH at high roughness**.

## The Real Issue

**I2EM's single transition factor doesn't work for AIEM's HH!**

The problem is:
1. I2EM uses **same Tf for both V and H** (0.35)
2. AIEM needs **different Tf for V and H** (0.81 vs 0.93)
3. HH specifically needs a **larger transition** than VV
4. Using I2EM's small Tf = 0.35 for HH causes failure at high roughness

## Why VV Works But HH Doesn't

**VV**: 
- AIEM original: Tfv = 0.81
- I2EM: Tf = 0.35
- Difference: 0.46 (moderate change)
- Result: VV improves! ✓

**HH**:
- AIEM original: Tfh = 0.93
- I2EM: Tf = 0.35
- Difference: 0.58 (huge change!)
- Result: HH fails at high roughness! ❌

The I2EM transition is **too aggressive a change for HH**.

## Solutions

### Option 1: Use I2EM Transition Only for VV

```python
if self._use_i2em_transition:
    Tfv, _ = compute_i2em_transition_function(...)
    _, Tfh = compute_transition_function(...)  # Keep AIEM for HH
```

**Pros**: VV improves, HH stays stable
**Cons**: Inconsistent approach

### Option 2: Scale I2EM Transition for HH

```python
if self._use_i2em_transition:
    Tfv, Tfh_i2em = compute_i2em_transition_function(...)
    # Scale up Tfh to match AIEM's ratio
    Tfh = Tfh_i2em * (0.93 / 0.35)  # Scale by AIEM's HH preference
```

**Pros**: Maintains I2EM philosophy but adapts for HH
**Cons**: Arbitrary scaling factor

### Option 3: Hybrid Transition (RECOMMENDED)

Use I2EM transition for VV, AIEM transition for HH:

```python
if self._use_i2em_transition:
    # Use I2EM for VV (works well)
    Tfv, _ = compute_i2em_transition_function(eps_r, theta_i, ks, cs, spectra, n_terms)
    
    # Use AIEM for HH (needed for high roughness)
    _, Tfh = compute_transition_function(eps_r, theta_i, ks, cs, spectra, n_terms)
```

This gives:
- ✅ VV: I2EM transition (RMSE ≈ 1.5-2.0 dB)
- ✅ HH: AIEM transition (RMSE ≈ 2-4 dB, stable at high roughness)
- ✅ Both work across full roughness range

### Option 4: Roughness-Dependent Transition

Blend I2EM and AIEM transitions based on roughness:

```python
if self._use_i2em_transition:
    Tfv_i2em, Tfh_i2em = compute_i2em_transition_function(...)
    Tfv_aiem, Tfh_aiem = compute_transition_function(...)
    
    # Blend based on ks
    alpha = min(1.0, ks / 0.8)  # 0 at low ks, 1 at high ks
    Tfv = (1 - alpha) * Tfv_i2em + alpha * Tfv_aiem
    Tfh = (1 - alpha) * Tfh_i2em + alpha * Tfh_aiem
```

**Pros**: Smooth transition, works at all roughness levels
**Cons**: More complex

## Recommendation

**Use Option 3: Hybrid Transition**

This is the simplest and most effective solution:
- VV gets I2EM's better transition (fixes the +5 dB bias)
- HH keeps AIEM's transition (stable at high roughness)
- No complex blending or scaling needed

## Expected Results with Hybrid Approach

| Band | Pol | Current (I2EM for both) | Hybrid (I2EM for VV, AIEM for HH) |
|------|-----|------------------------|-----------------------------------|
| 4.0 | VV | 2.11 dB ✓ | 2.11 dB ✓ (no change) |
| 4.0 | HH | 5.81 dB ❌ | **2-3 dB ✓** (much better) |
| 15.0 | VV | 1.58 dB ✓ | 1.58 dB ✓ (no change) |
| 15.0 | HH | 8.50 dB ❌ | **2-3 dB ✓** (much better) |

The hybrid approach should give:
- VV: RMSE ≈ 1.5-2.0 dB (I2EM quality)
- HH: RMSE ≈ 2-4 dB (AIEM quality, stable)
- Both work across full roughness range

## Why This Makes Physical Sense

**VV and HH have different scattering mechanisms:**

- **VV**: More sensitive to surface tilt, benefits from gentler transition
- **HH**: More sensitive to surface roughness, needs stronger transition at high ks

I2EM's single transition factor is a simplification that works for I2EM's simpler formulation, but AIEM's more complete physics requires different transitions for V and H polarizations.

**The hybrid approach respects this physical difference!**
