# IEM Family Refactor Plan

## Objective
Identify common building blocks across the MATLAB I2EM/IEM implementations and outline
how to encapsulate them into reusable components inside `mwrtms`, making it easy to add or
extend IEM variants (classic IEM, AIEM, I2EM, multiscale).

## Reference Sources
- `.temp/MATLAB/mrs.eecs.umich.edu/I2EM_Backscatter_model.m`
- `.temp/MATLAB/mrs.eecs.umich.edu/I2EM_Bistat_model.m`
- `.temp/MATLAB/mrs.eecs.umich.edu/IEMX_model.m`
- `.temp/MATLAB/mrs.eecs.umich.edu/Multiscale_I2EM_Backscatter.m`

## Core Building Blocks Observed
1. **Geometry and Fresnel Terms**
   - `theta`, `theta_s`, `phi_s` conversions.
   - Fresnel reflection coefficients `Rvi`, `Rhi`, plus average versions `Rav`, `Rah` for non-backscatter.
   - Shared between I2EM_Bistat and Multiscale variants.

2. **Surface Statistics**
   - Roughness parameters `ks = k * sigma`, `kl = k * L`.
   - Surface spectra (`roughness_spectrum`, `spectrm1/spectrm2`).
   - RMS slope `rss` derived from correlation function type (exponential, Gaussian, power-law).

3. **Shadowing Functions**
   - Smith (1967) shadowing (`ShdwS`, `Shdw`) used for both co-pol and cross-pol corrections.

4. **Series Expansion**
   - Number of spectral terms `Ts` selected via error tolerance.
   - Summations over `n` (and `m` for cross-pol) using factorial denominators and spectral weights.

5. **Cross-Polar Component (IEMX_model)**
   - Double integral evaluated via `dblquad` over `r, phi`.
   - Reuses same `ks`, `kl`, reflection coeffs but with different combination (`rvh` etc.).

## Proposed Abstractions
1. **`IEMKernel` Base Class** (Python)
   - Responsibilities:
     - Accept `surface`, `wave`, `geometry`, `dielectric`.
     - Provide cached Fresnel terms and slopes.
     - Provide helpers for spectral weight generation (exponential/Gaussian/power-law).
     - Provide Smith shadowing utility.
   - Methods to override: `compute_copol`, `compute_crosspol`.

2. **`IEMCorrelation` utilities**
   - Encapsulate correlation-specific calculations: `rss`, `roughness_spectrum`, `spectrm1/2`.
   - Likely implemented as strategy classes or functions reused by both base and cross-pol modules.

3. **`IEMSeries` helper**
   - Manage spectral order selection (auto vs fixed) and factorial sums.

4. **`IEMCrossPolar` integral module**
   - Translate `xpol_integralfunc`, `spectrm1`, `spectrm2` into vectorized NumPy/SciPy code.
   - Should support auto spectral component selection like MATLAB `auto=1`.

5. **`IEMResult` integration**
   - Standardize outputs via `ScatteringResult`.

## Implementation Steps
1. **Port Common Utilities**
   - Implement correlation helper functions matching MATLAB `roughness_spectrum`, `spectrm1/2`, etc.
   - Validate outputs against MATLAB by replicating small test cases.

2. **Create `IEMBase` class**
   - Provide `compute_co_pol()` returning HH/VV using shared formulas.
   - Provide hooking mechanism for AIEM-specific modifications (e.g., complementary term, slope probability adjustments).

3. **Implement Cross-Pol Module**
   - Port `IEMX_model` to Python, linking to base class utilities.
   - Ensure integration uses SciPy or fallback numeric integration (consider performance).

4. **Refactor Existing Models**
   - Update `SPMModel` (already compatible) as baseline data check.
   - Create new `I2EMModel` using base class.
   - Later: adapt `AIEMModel` to new base (reuse complementary calculations).

5. **Testing Strategy**
   - Generate fixtures from MATLAB scripts (e.g., call MATLAB code via saved CSV).
   - Add unit tests verifying co-pol and cross-pol results within tolerance (~0.1 dB).

6. **Documentation**
   - Describe new class hierarchy and extension points.
   - Provide example usage comparing I2EM vs AIEM.

## Progress
- ✅ Ported correlation helper functions into `mwrtms.scattering.iem.correlation`.
- ✅ Added `IEMBase` skeleton with shared single-scale preprocessing.
- ✅ Implemented `I2EMModel` (co-polar components) backed by `IEMBase`.

## Immediate Next Action
- Extend cross-polar computations (IEMX) and integrate validation datasets.
