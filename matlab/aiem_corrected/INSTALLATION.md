# Installation and Usage Guide

## Quick Start

### 1. Add to MATLAB Path

```matlab
% Add the aiem_corrected directory to your MATLAB path
addpath('/path/to/mwrtms/matlab/aiem_corrected');
```

Or use the MATLAB GUI: Home → Set Path → Add Folder

### 2. Basic Usage

```matlab
% Simple example
[VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.5, 15.0, 1.5, 2);
fprintf('VV = %.2f dB\n', VV);
```

### 3. Run Tests

```matlab
% Run the test suite
test_aiem_corrected
```

---

## Detailed Installation

### Requirements

- MATLAB R2016b or later (tested on R2020a+)
- No additional toolboxes required
- All functions use base MATLAB only

### File Structure

```
matlab/aiem_corrected/
├── AIEM_corrected.m              % Main function
├── aiem_single_scattering.m      % Single scattering core
├── fresnel_coefficients.m        % Fresnel coefficients (bug-fixed)
├── transition_function.m         % Transition function (bug-fixed)
├── kirchhoff_coefficients.m      % Kirchhoff field coefficients
├── complementary_coefficients.m  % Complementary field (bug-fixed)
├── roughness_spectrum.m          % Roughness spectrum (bug-fixed)
├── test_aiem_corrected.m         % Test suite
├── compare_with_nmm3d.m          % NMM3D comparison
├── README.md                     % Main documentation
├── BUG_FIXES.md                  % Detailed bug documentation
└── INSTALLATION.md               % This file
```

---

## Usage Examples

### Example 1: Basic Monostatic Backscatter

```matlab
% Parameters
theta_i = 40;      % Incident angle (degrees)
theta_s = 40;      % Scattered angle (degrees)
phi_s = 180;       % Backscatter azimuth (degrees)
kl = 5.0;          % Normalized correlation length
ks = 0.5;          % Normalized RMS height
err = 15.0;        % Real part of permittivity
eri = 1.5;         % Imaginary part of permittivity
itype = 2;         % Exponential correlation

% Compute backscatter
[VV, HH, HV, VH] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, err, eri, itype);

% Display results
fprintf('Backscatter coefficients:\n');
fprintf('  VV = %.2f dB\n', VV);
fprintf('  HH = %.2f dB\n', HH);
fprintf('  HV = %.2f dB\n', HV);
fprintf('  VH = %.2f dB\n', VH);
```

### Example 2: Angular Scan

```matlab
% Parameters
angles = 20:5:60;
kl = 5.0;
ks = 0.4;
err = 12.0;
eri = 1.0;
itype = 2;

% Compute for each angle
VV_array = zeros(size(angles));
HH_array = zeros(size(angles));

for i = 1:length(angles)
    theta = angles(i);
    [VV, HH, ~, ~] = AIEM_corrected(theta, theta, 180, kl, ks, err, eri, itype);
    VV_array(i) = VV;
    HH_array(i) = HH;
end

% Plot results
figure;
plot(angles, VV_array, 'b-o', 'LineWidth', 2);
hold on;
plot(angles, HH_array, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Incidence Angle (degrees)');
ylabel('Backscatter (dB)');
legend('VV', 'HH');
title('Angular Dependence');
```

### Example 3: Roughness Scan

```matlab
% Parameters
theta = 40;
kl = 5.0;
ks_values = 0.1:0.1:1.0;
err = 15.0;
eri = 1.5;
itype = 2;

% Compute for each roughness
VV_array = zeros(size(ks_values));
HH_array = zeros(size(ks_values));

for i = 1:length(ks_values)
    ks = ks_values(i);
    [VV, HH, ~, ~] = AIEM_corrected(theta, theta, 180, kl, ks, err, eri, itype);
    VV_array(i) = VV;
    HH_array(i) = HH;
end

% Plot results
figure;
plot(ks_values, VV_array, 'b-o', 'LineWidth', 2);
hold on;
plot(ks_values, HH_array, 'r-s', 'LineWidth', 2);
grid on;
xlabel('Normalized RMS Height (ks)');
ylabel('Backscatter (dB)');
legend('VV', 'HH');
title('Roughness Dependence');
```

### Example 4: Different Correlation Types

```matlab
% Parameters
theta = 40;
kl = 5.0;
ks = 0.5;
err = 15.0;
eri = 1.5;

% Correlation types
corr_names = {'Gaussian', 'Exponential', '1.5-power'};

fprintf('Correlation Type    VV (dB)    HH (dB)\n');
fprintf('----------------------------------------\n');

for itype = 1:3
    [VV, HH, ~, ~] = AIEM_corrected(theta, theta, 180, kl, ks, err, eri, itype);
    fprintf('%-15s  %7.2f    %7.2f\n', corr_names{itype}, VV, HH);
end
```

### Example 5: Physical Parameters to Normalized

```matlab
% Physical parameters
frequency_ghz = 5.405;        % Frequency (GHz)
rms_height_cm = 1.0;          % RMS height (cm)
corr_length_cm = 5.0;         % Correlation length (cm)
soil_moisture = 0.25;         % Volumetric moisture (0-1)

% Convert to normalized parameters
lambda_m = 0.3 / frequency_ghz;  % Wavelength (m)
k = 2 * pi / lambda_m;           % Wavenumber (rad/m)
ks = k * (rms_height_cm / 100);  % Normalized RMS height
kl = k * (corr_length_cm / 100); % Normalized correlation length

% Estimate permittivity (simple model)
eps_r_real = 3.0 + 20.0 * soil_moisture;
eps_r_imag = 0.5 + 2.0 * soil_moisture;

% Compute backscatter
theta = 40;
itype = 2;
[VV, HH, HV, VH] = AIEM_corrected(theta, theta, 180, kl, ks, eps_r_real, eps_r_imag, itype);

fprintf('Physical parameters:\n');
fprintf('  Frequency: %.3f GHz\n', frequency_ghz);
fprintf('  RMS height: %.2f cm\n', rms_height_cm);
fprintf('  Correlation length: %.2f cm\n', corr_length_cm);
fprintf('  Soil moisture: %.2f\n\n', soil_moisture);

fprintf('Normalized parameters:\n');
fprintf('  ks = %.4f\n', ks);
fprintf('  kl = %.4f\n\n', kl);

fprintf('Results:\n');
fprintf('  VV = %.2f dB\n', VV);
fprintf('  HH = %.2f dB\n', HH);
```

---

## Validation

### Run Test Suite

```matlab
% Run all tests
test_aiem_corrected
```

Expected output:
```
✓ All bug fixes applied
✓ Fresnel branch correction working
✓ Normal incidence constants correct
✓ Spectrum scaling laws verified
✓ Monostatic reciprocity satisfied
```

### Compare with NMM3D

```matlab
% Requires NMM3D reference data
compare_with_nmm3d
```

Expected performance:
- VV: RMSE ≈ 3 dB, Bias ≈ +3 dB
- HH: RMSE ≈ 5 dB, Bias ≈ +5 dB

---

## Troubleshooting

### Issue: "Function not found"

**Solution:** Add directory to MATLAB path:
```matlab
addpath('/path/to/matlab/aiem_corrected');
savepath;  % Save for future sessions
```

### Issue: "NMM3D data file not found"

**Solution:** Update path in `compare_with_nmm3d.m`:
```matlab
nmm3d_file = 'path/to/NMM3D_LUT_NRCS_40degree.dat';
```

### Issue: Results differ from Python

**Expected:** Small numerical differences (<0.01 dB) due to:
- Floating-point precision
- Different implementations of special functions
- Convergence criteria

**Not expected:** Large differences (>0.1 dB) - check input parameters

### Issue: Very large or -Inf values

**Causes:**
- HV/VH with single scattering only (expected for smooth surfaces)
- Invalid input parameters (ks < 0, kl < 0, etc.)
- Extreme angles (>70 degrees)

**Solution:** Check input validity and consider multiple scattering for cross-pol

---

## Performance Tips

### For Batch Processing

```matlab
% Pre-allocate arrays
n_configs = 1000;
results = zeros(n_configs, 4);  % [VV, HH, HV, VH]

% Vectorize parameter generation
ks_array = linspace(0.1, 1.0, n_configs);

% Compute in loop (cannot vectorize AIEM itself)
for i = 1:n_configs
    [VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, ks_array(i), 15.0, 1.5, 2);
    results(i, :) = [VV, HH, HV, VH];
end
```

### For Speed

- Use `itype = 2` (Exponential) - fastest
- Avoid `itype = 3` (1.5-power) - slower due to more complex spectrum
- Reduce `ks` and `kl` if possible - fewer spectral terms needed

---

## Known Limitations

1. **Systematic bias:** +3-5 dB vs NMM3D for co-pol
   - Root cause: Legacy transition function
   - Solution: Implement new S_p/S_p^(0) method

2. **Cross-pol:** Single scattering only
   - HV/VH often -Inf or very small
   - Solution: Add multiple scattering module

3. **Large angles:** Accuracy decreases above 60°
   - Consider adding shadowing function

4. **Very rough surfaces:** ks > 3 may be inaccurate
   - AIEM validity range: ks < 3

---

## Support

For issues or questions:
1. Check README.md and BUG_FIXES.md
2. Review test_aiem_corrected.m for examples
3. See Python implementation for reference
4. Consult original papers (Chen et al. 2003, Wu & Fung 1992)

---

## License

Same as parent project.

---

**Version:** 1.0 (2024)  
**Status:** ✅ All bugs fixed | ⚠️ Legacy transition function
