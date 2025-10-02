# mwRTMs

mwRTMs (Microwave Radiative Transfer Models) is a Python toolbox that groups
active and passive microwave radiative transfer models under a single
object-oriented architecture. The codebase is structured around Python's four
OOP pillars—encapsulation, inheritance, polymorphism, and abstraction—so new
models can plug into shared interfaces without sacrificing clarity or
extensibility.

---

## Features

- **Unified OOP architecture**: Core types such as
  `ElectromagneticWave`, `ScatteringGeometry`, `SurfaceRoughness`, and
  `Medium` enforce unit conversions, caching, and validation so solver
  implementations focus on the physics. Subclasses inherit the base
  `SurfaceScattering` and `VolumeScattering` templates to encourage
  encapsulation and polymorphism.
- **Surface RTM portfolio**: AIEM, PRISM1, Dubois95, SMART, SPM3D, and I2EM are
  provided via dedicated surface-scattering classes that adhere to the common
  mwRTMs interfaces. (Solver implementations are currently disabled while the
  codebase is refactored.)
- **Vegetation RTMs**: The single-scattering radiative transfer (SSRT) canopy
  wrapper remains available in the volume-scattering hierarchy, with the
  underlying solver temporarily disabled during the same refactor.
- **Reference utilities**: Helper routines (Dobson85 dielectric model, Fresnel
  coefficients, dB/power conversions) stay accessible through the mwRTMs
  utility modules for reuse in custom workflows.

---

## Requirements

- Python 3.9+
- NumPy, SciPy, Matplotlib, Pandas (see `requirements.txt` for the full list)
- A C/C++ toolchain to build the optional I2EM Cython extension (`build-essential`
  on Debian/Ubuntu, Xcode Command Line Tools on macOS, or Visual Studio Build
  Tools on Windows).

---

## Installation

```bash
# clone the repository
git clone https://github.com/<your-org>/mwRTMs.git
cd mwRTMs

# (optional) create a virtual environment
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# install the package in editable mode
pip install -r requirements.txt
pip install --upgrade pip setuptools Cython numpy
pip install -e .
```

You can also install directly from GitHub using pip:

```bash
pip install git+https://github.com/Morteza-Khazaei/mwRTMs.git
```

Thanks to `pyproject.toml`, pip will pull in the build requirements
(`setuptools`, `wheel`, `Cython`, and `numpy>=1.21`) before compiling the I2EM
Cython extension. If the compiler toolchain is unavailable pip automatically
falls back to the pure-Python I2EM backend.

---

## Quick Start

```python
from mwrtms import (
    ElectromagneticWave,
    IsotropicMedium,
    PolarizationState,
    ScatteringGeometry,
    SurfaceRoughness,
    SSRTModel,
)

wave = ElectromagneticWave(frequency_ghz=5.405)
geometry = ScatteringGeometry(theta_i_deg=40.0, theta_s_deg=40.0, phi_s_deg=180.0)
roughness = SurfaceRoughness(rms_height=0.02, correlation_length=0.1)

canopy = IsotropicMedium(permittivity=12.0 + 3.0j)
soil = IsotropicMedium(permittivity=5.0 + 1.0j)

ssrt = SSRTModel(
    wave,
    geometry,
    roughness,
    single_scatter_albedo=0.08,
    extinction_coefficient=0.5,
    canopy_thickness=0.3,
    surface_model="I2EM",
    canopy_model="Diff",
)

sigma_vv = ssrt.backscatter(canopy, soil, polarization=PolarizationState.VV)
print(f"Total VV backscatter: {sigma_vv:.3f} [linear]")
```

---

## Example Notebooks

Located under `examples/`:

| Notebook | Description |
| --- | --- |
| `test_spm.ipynb` | SPM3D vs NMM3D LUT comparisons grouped by ℓ/σ. |
| `test_prism1.ipynb` | PRISM1 benchmark against NMM3D. |
| `test_smart.ipynb` | SMART validation and soil sweep. |
| `test_i2em.ipynb` | I2EM (bistatic) interpolation example. |
| `test_dubois95.ipynb` | Dubois95 LUT comparison. |
| `test_ssrt.ipynb` | Soil–vegetation integration demo with incidence/extinction sweeps. |

Launch a notebook with:

```bash
jupyter notebook examples/test_ssrt.ipynb
```

---

## Testing & Validation

1. Run the unit-style scripts in `examples/` after modifying a soil model.
2. Use the LUT notebooks to confirm backscatter errors remain within acceptable
   tolerances (≈1 dB for VV/HH).
3. For vegetation scenarios, rerun the `test_ssrt.ipynb` canopy sweeps to gauge
   sensitivity to new parameter ranges.

---

## Contributing

1. Fork the repository and create a feature branch.
2. Add or update notebooks/tests for any new model or parameterisation.
3. Run linting/tests and ensure notebooks execute cleanly.
4. Submit a pull request summarising the change and validation steps.

---

## License

This project is released under the MIT License. See `LICENSE` for details.

---

## Citing

If you use mwRTMs in academic work, please cite the repository and the relevant
model references (AIEM, PRISM1, Dubois95, SMART, SPM3D, I2EM, SSRT) as
appropriate.
