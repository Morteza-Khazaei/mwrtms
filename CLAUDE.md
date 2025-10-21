# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**mwRTMs** (Microwave Radiative Transfer Models) is a Python library for active microwave backscatter modeling focused on soil and vegetation. It implements advanced electromagnetic scattering models with clean OOP abstractions and Numba acceleration for performance.

## Development Commands

### Installation
```bash
pip install -e .
```

### Running Tests
Tests are located in `tests/` directory. Run individual test scripts directly:
```bash
# Validation tests against NMM3D reference data
PYTHONPATH=src python3 tests/aiem_nmm3d_test.py
PYTHONPATH=src python3 tests/i2em_nmm3d_test.py
PYTHONPATH=src python3 tests/ka_nmm3d_test.py

# Performance benchmarks
PYTHONPATH=src python3 tests/benchmark_aiem_numba.py

# Unit tests with pytest (if configured)
PYTHONPATH=src python3 -m pytest tests/unit/
```

### Running Examples
Interactive examples are in Jupyter notebooks:
```bash
# Navigate to examples directory
cd examples/

# Run notebooks
jupyter notebook test_aiem.ipynb
jupyter notebook test_ka.ipynb
jupyter notebook test_spm.ipynb
```

## Architecture Overview

### Core Abstractions (src/mwrtms/core/)

The package is built on fundamental electromagnetic abstractions:
- **ElectromagneticWave**: Wave properties (frequency, wavelength, wavenumber)
- **ScatteringGeometry**: Incident/scattered angles (monostatic and bistatic)
- **PolarizationState**: HH, VV, HV, VH polarization channels
- **RadarConfiguration**: Observation modes (monostatic/bistatic)

### Medium Representations (src/mwrtms/medium/)

Physical media are represented through:
- **Medium**: Base class for dielectric media
- **SoilMedium/MironovSoilMedium**: Soil dielectric properties
- **VegetationMedium**: Vegetation canopy properties
- **Surface**: Surface roughness (height statistics, correlation)
- **FresnelCoefficients**: Reflection/transmission at interfaces

### Scattering Models (src/mwrtms/scattering/)

Implements multiple surface scattering models as separate classes:

**Surface Models:**
- **AIEMModel** (Advanced IEM): Primary model with multiple scattering support
  - Located in `scattering/surface/iem/aiem.py`
  - Includes Kirchhoff, complementary, and multiple scattering terms
  - Supports transition functions for smooth I2EM convergence
  - Has Numba-accelerated backend in `aiem_numba_backend.py`
  - Multiple scattering implementation in `multiple_scattering.py`

- **I2EMModel** (Improved IEM): Legacy implementation
  - Located in `scattering/surface/iem/i2em.py`

- **SPMModel** (Small Perturbation Method): For very smooth surfaces
  - Located in `scattering/surface/spm.py`

- **KirchhoffApproximation**: For very rough surfaces
  - Located in `scattering/surface/ka.py`

**Volume Models:**
- **SSRTModel**: Simple first-order volume scattering
  - Located in `scattering/volume/`

### Factory Pattern (src/mwrtms/factory/)

Models register themselves using `@register_model` decorator:
```python
@register_model("aiem")
class AIEMModel(IEMBase):
    ...
```

Allows creation via:
```python
model = ScatteringModelFactory.create("aiem", wave, geometry, surface)
```

### Facade API (src/mwrtms/facade/)

Simplified high-level API through `mwRTMs` class provides user-friendly interface for common operations.

## IEM Implementation Details

### AIEM Multiple Scattering Architecture

The AIEM model implements second-order multiple scattering for cross-polarization (HV/VH):

**Key Components:**
1. **Propagators** (Fp, Fm, Gp, Gm): Upward/downward electromagnetic fields
2. **Kirchhoff-Complementary Terms** (K1, K2, K3): From Yang et al. (2017)
3. **Complementary Terms** (gc1-gc14): 14 field interaction terms
4. **2D Spectral Integration**: Trapezoidal quadrature over (U,V) spectral domain

**Numba Acceleration:**
- Multiple scattering uses Numba JIT compilation for ~100x speedup
- Fallback to pure NumPy if Numba not available
- First call is slow (~2s) due to JIT compilation
- Subsequent calls are fast (~0.17s with standard settings)

**Guardrails System:**
Located in `scattering/surface/iem/guardrails.py`:
- Input validation (wavenumbers, angles, roughness)
- Fresnel coefficient bounds checking
- Energy conservation verification
- Cross-polarization Kirchhoff term validation
- Multiple scattering balance checks

### Critical Implementation Notes

1. **Spectrum Computation**:
   - nth-order spectrum W^(n) used in double summations
   - Exponential correlation: W^(n)(k) = 4πℓ²n / [n² + (ℓk)²]^(3/2)
   - Must include σ² (RMS height squared) in normalization

2. **Fresnel Coefficients**:
   - Cross-pol uses R ≈ (Rh + Rv)/2 approximation
   - Complex wavenumber in lower medium: kt_z has specific sign conventions
   - Transition functions smooth convergence to I2EM at small roughness

3. **Multiple Scattering Integration**:
   - Domain: typically [-10/kℓ, 10/kℓ] for adequate coverage
   - Quadrature points: 65 (fast), 129 (standard), 257 (high-res)
   - Radiation condition masking applied to avoid evanescent waves
   - Prefactors: k²/(8π) for copol terms, k²/(64π) for cross-pol

4. **Known Limitations**:
   - HV has ~31 dB systematic bias vs NMM3D (use for relative comparisons)
   - VV accuracy: <3 dB RMSE
   - HH accuracy: <5 dB RMSE

## Testing Strategy

### Validation Approach

Tests compare against **NMM3D** (Numerical Method of Moments in 3D) reference data:
- NMM3D provides ground truth for surface scattering
- Test data stored in `data/NMM3D_LUT_NRCS_*.dat` files
- 162 test cases covering various roughness conditions
- Metrics: RMSE, MAE, bias, Pearson correlation

### Test Organization

```
tests/
├── *_nmm3d_test.py      # Validation against NMM3D
├── benchmark_*.py       # Performance benchmarks
├── debug_*.py           # Diagnostic scripts
├── diagnose_*.py        # Problem analysis
├── unit/                # Unit tests
│   └── surface/         # Surface scattering tests
└── conftest.py          # pytest configuration
```

### Running Validation Tests

```bash
# AIEM validation with multiple scattering
PYTHONPATH=src python3 tests/aiem_nmm3d_test.py --add-multiple

# I2EM validation
PYTHONPATH=src python3 tests/i2em_nmm3d_test.py

# Kirchhoff approximation validation
PYTHONPATH=src python3 tests/ka_nmm3d_test.py
```

## Documentation

Extensive technical documentation in `docs/`:

**Quick References:**
- `QUICK_START.md` - 5-minute getting started guide
- `FINAL_SUMMARY.md` - Complete implementation summary

**Technical Specifications:**
- `AIEM_MULTIPLE_SCATTERING_INSTRUCTION.md` - Mathematical blueprint
- `aiem_multiple_scattering_spec.md` - Condensed specification
- `NUMBA_ACCELERATION_GUIDE.md` - Performance optimization

**Implementation Notes:**
- `GUARDRAILS_*.md` - Validation system documentation
- `KA_*.md` - Kirchhoff approximation notes
- Various `AIEM_MS_*.md` - Development history and bug fixes

## Key Dependencies

- **numpy**: Array operations and numerical computing
- **scipy**: Integration, special functions
- **numba** (optional): JIT compilation for 100x speedup
- **matplotlib**: Plotting (examples)
- **pandas**: Data handling (examples)

Python ≥3.10 required.

## Common Patterns

### Creating and Running a Model

```python
from mwrtms import AIEMModel, ElectromagneticWave, ScatteringGeometry
from mwrtms import build_surface_from_statistics, HomogeneousMedium

# Setup geometry and wave
wave = ElectromagneticWave(5.4e9)  # 5.4 GHz
geometry = ScatteringGeometry(theta_i_deg=40.0)

# Define surface
surface = build_surface_from_statistics(
    rms_height=0.01,           # 1 cm
    correlation_length=0.05,   # 5 cm
    correlation_type="exponential"
)

# Create model
model = AIEMModel(
    wave, geometry, surface,
    include_multiple_scattering=True  # Enable for HV/VH
)

# Define media
air = HomogeneousMedium(1.0 + 0.0j)
soil = HomogeneousMedium(12.0 + 1.8j)

# Run computation
result = model.run(air, soil)

# Access results
print(f"VV: {result.vv_db:.2f} dB")
print(f"HH: {result.hh_db:.2f} dB")
print(f"HV: {result.hv_db:.2f} dB")
```

### Model Configuration Options

```python
# Fast mode (testing)
model = AIEMModel(
    wave, geometry, surface,
    include_multiple_scattering=True,
    ms_quadrature_points=65,
    ms_spectral_terms=6
)

# Standard mode (recommended)
model = AIEMModel(
    wave, geometry, surface,
    include_multiple_scattering=True,
    ms_quadrature_points=129,
    ms_spectral_terms=8
)

# High-resolution mode (publication)
model = AIEMModel(
    wave, geometry, surface,
    include_multiple_scattering=True,
    ms_quadrature_points=257,
    ms_spectral_terms=10
)
```

## Git Workflow Notes

Based on current git status, the repository shows:
- Active development on AIEM multiple scattering implementation
- Several staged/unstaged changes to IEM models
- New guardrails system being integrated
- Extensive documentation being updated

When committing changes:
- Focus commits on single logical changes
- Reference relevant documentation files
- Include test validation if modifying scattering models
- Update corresponding docs in `docs/` directory
