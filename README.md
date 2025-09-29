# SSRT

Single Scattering Radiative Transfer (SSRT) is a Python toolbox for simulating radar backscatter from vegetated surfaces. It couples first-order surface scattering models with a layered canopy representation so you can explore how soil dielectric, roughness, and vegetation structure jointly control the microwave response.

---

## Features

- **Soil backscatter models**: AIEM, PRISM1, Dubois95, SMART, SPM3D, and I2EM are exposed through a common interface and can be mixed with vegetation layers via `S2RTR`.
- **Canopy integration**: The `S2RTR` core supports diffuse and specular upper-boundary formulations with configurable albedo, extinction coefficient, thickness, and canopy permittivity.
- **Reference benchmarking**: Ready-to-run Jupyter notebooks validate each soil model against the 40° NMM3D LUT and illustrate SSRT soil–vegetation integration scenarios.
- **Configurable incidence geometry**: Control incidence and scattering angles, azimuth, and correlation statistics to match experimental or mission conditions.

---

## Requirements

- Python 3.9+
- NumPy, SciPy, Matplotlib, Pandas (see `requirements.txt` for the full list)
- A C/C++ toolchain when installing from source so the optional Cython extension can compile (e.g. `build-essential` on Debian/Ubuntu, Xcode Command Line Tools on macOS, or Visual Studio Build Tools on Windows).

---

## Installation

```bash
# clone the repository
 git clone https://github.com/<your-org>/SSRT.git
 cd SSRT

# (optional) create a virtual environment
 python -m venv .venv
 source .venv/bin/activate  # Windows: .venv\Scripts\activate

# install the package in editable mode
 pip install -r requirements.txt
 pip install -e .
```

You can also install directly from GitHub using pip:

```bash
pip install git+https://github.com/Morteza-Khazaei/SSRT.git
```

Thanks to `pyproject.toml`, pip will pull in the build requirements (`setuptools`, `wheel`, `Cython`, and `numpy>=1.21`) before compiling the I2EM Cython extension. Ensure that the system compiler toolchain mentioned above is available; otherwise pip will fall back to the pure-Python implementation of I2EM.

---

## Quick Start

```python
import numpy as np
from ssrt import S2RTR

soil_eps = np.array([5.0 + 2.0j])
canopy_eps = np.array([12.0 + 3.0j])

rt = S2RTR(
    frq_GHz=5.405,
    theta_i=40.0,
    theta_s=40.0,
    phi_i=0.0,
    phi_s=180.0,
    s=0.03,          # soil RMS height (m)
    cl=0.12,         # correlation length (m)
    eps2=canopy_eps,
    eps3=soil_eps,
    a=0.08,
    kappa_e=0.7,     # extinction coefficient (Np/m)
    d=0.35,          # canopy thickness (m)
    acftype='exp',   # soil autocorrelation (exp, gauss, pow)
    RT_models={'RT_s': 'SPM3D', 'RT_c': 'Spec'}
)

sig_ground, sig_soil, sig_total = rt.calc_sigma(todB=True)
print('Ground VV (dB):', sig_soil['vv'])
print('Total VV (dB):', sig_total['vv'])
```

---

## Example Notebooks

Located under `test/`:

| Notebook | Description |
| --- | --- |
| `test_spm.ipynb` | SPM3D vs NMM3D LUT comparisons grouped by ℓ/σ. |
| `test_prism1.ipynb` | PRISM1 benchmark against NMM3D. |
| `test_smart.ipynb` | SMART validation and soil sweep. |
| `test_i2em.ipynb` | I2EM (bistatic) interpolation example. |
| `test_dubois95.ipynb` | Dubois95 LUT comparison. |
| `test_ssrt.ipynb` | Soil–vegetation integration demo with incidence/extinction sweeps. |

To launch:

```bash
jupyter notebook test/test_ssrt.ipynb
```

---

## Testing & Validation

1. Run unit-style scripts or notebooks each time you modify a soil model.
2. Use the LUT notebooks to confirm backscatter errors remain within acceptable tolerances (≈1 dB for VV/HH).
3. For vegetation scenarios, rerun the `test_ssrt.ipynb` canopy sweeps to gauge sensitivity to new parameter ranges.

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

If you use SSRT in academic work, please cite the repository and the relevant model references (AIEM, PRISM1, Dubois95, SMART, SPM3D, I2EM) as appropriate.
