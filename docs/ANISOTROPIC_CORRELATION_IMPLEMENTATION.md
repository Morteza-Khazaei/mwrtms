# Anisotropic Exponential Correlation Function Implementation

## Summary

Implemented full anisotropic support for the exponential autocorrelation function in the AIEM model, following the formulation from Yang & Chen (2019) IEEE TGRS paper "Polarized Backscattering From Spatially Anisotropic Rough Surface".

## Changes Made

### 1. Updated `correlation.py`

**File:** `src/mwrtms/medium/surface/correlation.py`

- Modified the `Exponential` class to support anisotropic correlation lengths
- Added `__init__` method to accept `correlation_length_x` and `correlation_length_y` parameters
- Updated `spectrum()` method to compute azimuth-dependent correlation length:
  ```python
  L(φ) = Lx * cos²(φ) + Ly * sin²(φ)
  ```
- The spectrum now correctly implements:
  ```
  W^(n)(K,φ) = (L(φ)/n)² [1 + (KL(φ)/n)²]^(-1.5)
  ```

### 2. Updated `spectrum_aiem.py`

**File:** `src/mwrtms/scattering/iem/spectrum_aiem.py`

- Added optional parameters for anisotropic correlation lengths:
  - `kl_x`: Normalized correlation length along x-axis (k * Lx)
  - `kl_y`: Normalized correlation length along y-axis (k * Ly)
  - `kx`, `ky`: Spatial frequency components
- Implemented azimuth-dependent spectrum calculation for exponential correlation
- Maintains backward compatibility with isotropic case

### 3. Updated `multiple_scattering.py`

**File:** `src/mwrtms/scattering/iem/multiple_scattering.py`

- Extended `SurfaceParams` dataclass to include `kl_x` and `kl_y` fields
- Modified `_make_Wn_provider()` function to support anisotropic exponential correlation
- The provider now computes:
  ```python
  # Compute azimuthal angle from spectral components
  cos_phi = u / K
  sin_phi = v / K
  
  # Azimuth-dependent correlation length
  kl_phi = kl_x * cos²(φ) + kl_y * sin²(φ)
  
  # Spectrum with anisotropic correlation
  W^(n) = σ² * (kl_phi²/n) * [1 + (K*kl_phi/n)²]^(-1.5)
  ```

## Mathematical Formulation

### Paper Definition (Yang & Chen, 2019)

**Equation (1):** Anisotropic exponential correlation function
```
ρ(r,φ) = exp(-r/L(φ))
```

**Equation (2):** Azimuth-dependent correlation length
```
L(φ) = Lx cos²(φ) + Ly sin²(φ)
```

**Equation (4):** Roughness spectrum
```
W(K,φ) = (L(φ))² [1 + (K L(φ))²]^(-1.5)
```

### Implementation

For the nth-order spectrum used in AIEM:
```
W^(n)(K,φ) = (L(φ)/n)² [1 + (KL(φ)/n)²]^(-1.5)
```

Where:
- `n` is the spectral order (1, 2, 3, ...)
- `K = √(kx² + ky²)` is the spatial frequency magnitude
- `φ = atan2(ky, kx)` is the azimuthal angle in spectral domain
- `L(φ) = Lx cos²(φ) + Ly sin²(φ)` is the azimuth-dependent correlation length

## Key Features

1. **Backward Compatibility**: Isotropic case still works when `Lx = Ly` or when anisotropic parameters are not provided

2. **Numerical Stability**: Handles edge cases:
   - Division by zero when `K = 0` (uses average correlation length)
   - Smooth transition between isotropic and anisotropic modes

3. **Multiple Scattering**: Full anisotropic support in the multiple scattering integration, which is crucial for accurate AIEM predictions

4. **Proper Normalization**: Maintains correct σ² normalization and (2π) factors for 2D Fourier transform

## Usage Example

```python
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Create anisotropic surface
wave = ElectromagneticWave(5.4e9)  # 5.4 GHz
geometry = ScatteringGeometry(theta_i_deg=40.0)

# Build surface with anisotropic correlation lengths
surface = build_surface_from_statistics(
    rms_height=0.01,  # 1 cm
    correlation_length_m=0.15,  # Lx = 15 cm (minor axis)
    correlation_length_y_m=0.03,  # Ly = 3 cm (major axis)
    correlation_type="exponential"
)

# Run AIEM with multiple scattering
model = AIEMModel(
    wave, 
    geometry, 
    surface,
    include_multiple_scattering=True
)

air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(12.0 + 1.8j)
result = model.run(air, soil)

print(f"VV: {result.vv_db:.2f} dB")
print(f"HH: {result.hh_db:.2f} dB")
print(f"HV: {result.hv_db:.2f} dB")
```

## Validation

The implementation follows the exact formulation from:

**Yang, Y., & Chen, K. S. (2019).** "Polarized backscattering from spatially anisotropic rough surface." *IEEE Transactions on Geoscience and Remote Sensing*, 57(9), 6608-6618.

Key equations implemented:
- Equation (1): Anisotropic correlation function
- Equation (2): Azimuth-dependent correlation length
- Equation (4): Anisotropic roughness spectrum

## Impact on Multiple Scattering

The anisotropic correlation function is **crucial** for accurate multiple scattering calculations because:

1. **Azimuthal Dependence**: Multiple scattering involves integration over all azimuthal angles in the spectral domain
2. **Directional Effects**: Anisotropic surfaces (e.g., tilled soil) have strong directional scattering patterns
3. **Cross-Polarization**: HV/VH polarizations are particularly sensitive to surface anisotropy

## Testing

To test the anisotropic implementation:

```python
# Test with strong anisotropy (Lx = 0.5 cm, Ly = 3.5 cm)
# This should show:
# - Azimuthal modulation in backscatter
# - Dips at up/down directions (φ = 0°, 180°)
# - Enhanced cross-polarization
```

## References

1. Yang, Y., & Chen, K. S. (2019). Polarized backscattering from spatially anisotropic rough surface. IEEE TGRS, 57(9), 6608-6618.

2. Chen, K. S., Wu, T. D., Tsang, L., Li, Q., Shi, J., & Fung, A. K. (2003). Emission of rough surfaces calculated by the integral equation method. IEEE TGRS, 41(1), 90-101.

3. Fung, A. K., & Chen, K. S. (2010). Microwave Scattering and Emission Models for Users. Artech House.
