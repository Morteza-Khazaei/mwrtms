# AIEM Corrected - Complete Index

**Corrected MATLAB Implementation of AIEM with All Bug Fixes**

---

## 📁 File Organization

### Core Implementation
- **`AIEM_corrected.m`** - Main function (drop-in replacement for AIEM.m)
- **`aiem_single_scattering.m`** - Single scattering computation engine
- **`fresnel_coefficients.m`** - Fresnel reflection coefficients (✅ bug-fixed)
- **`transition_function.m`** - Transition function (✅ bug-fixed)
- **`kirchhoff_coefficients.m`** - Kirchhoff field coefficients
- **`complementary_coefficients.m`** - Complementary field coefficients (✅ bug-fixed)
- **`roughness_spectrum.m`** - Roughness spectrum computation (✅ bug-fixed)

### Testing & Validation
- **`test_aiem_corrected.m`** - Comprehensive test suite
- **`compare_with_nmm3d.m`** - Comparison with NMM3D reference data

### Documentation
- **`README.md`** - Main documentation and overview
- **`BUG_FIXES.md`** - Detailed bug fix documentation
- **`INSTALLATION.md`** - Installation and usage guide
- **`INDEX.md`** - This file

---

## 🎯 Quick Navigation

### For New Users
1. Start with **`README.md`** - Overview and quick start
2. Read **`INSTALLATION.md`** - Setup and examples
3. Run **`test_aiem_corrected.m`** - Verify installation

### For Developers
1. Review **`BUG_FIXES.md`** - Understand what was fixed
2. Study **`aiem_single_scattering.m`** - Core algorithm
3. Check **`compare_with_nmm3d.m`** - Validation approach

### For Researchers
1. See **`BUG_FIXES.md`** - Scientific justification
2. Review **`README.md`** - Performance vs NMM3D
3. Consult Python implementation for reference

---

## 📊 Bug Fixes Summary

| # | Bug | File | Status |
|---|-----|------|--------|
| 1 | Fresnel branch for lossy media | `fresnel_coefficients.m` | ✅ Fixed |
| 2 | Normal-incidence constants | `fresnel_coefficients.m`, `transition_function.m` | ✅ Fixed |
| 3 | Transition function polarization | `transition_function.m` | ✅ Fixed |
| 4 | 1.5-power spectrum | `roughness_spectrum.m` | ✅ Fixed |
| 5 | Complex magnitude checks | `complementary_coefficients.m` | ✅ Fixed |

**All bugs fixed!** ✅

---

## 🚀 Quick Start

```matlab
% Add to path
addpath('/path/to/matlab/aiem_corrected');

% Basic usage
[VV, HH, HV, VH] = AIEM_corrected(40, 40, 180, 5.0, 0.5, 15.0, 1.5, 2);
fprintf('VV = %.2f dB, HH = %.2f dB\n', VV, HH);

% Run tests
test_aiem_corrected
```

---

## 📈 Performance

### Current (with bug fixes)
- **VV:** RMSE = 2.93 dB, Bias = +2.77 dB
- **HH:** RMSE = 4.89 dB, Bias = +4.76 dB
- **Correlation:** >0.97 for both

### Known Limitation
- Systematic +3-5 dB bias due to legacy transition function
- **Root cause:** Using Wu & Fung (1992) method
- **Solution:** Implement new S_p/S_p^(0) method
- **Expected after fix:** RMSE < 1 dB

---

## 📚 Documentation Map

### Getting Started
1. **README.md** → Overview, features, quick start
2. **INSTALLATION.md** → Setup, examples, troubleshooting
3. **test_aiem_corrected.m** → Verify everything works

### Understanding the Fixes
1. **BUG_FIXES.md** → Detailed explanation of each bug
2. **Python implementation** → Reference implementation
3. **Bug report** → `../../docs/AIEM_MATLAB_BUG_REPORT.md`

### Validation
1. **compare_with_nmm3d.m** → Quantitative comparison
2. **test_aiem_corrected.m** → Sanity checks
3. **Root cause analysis** → `../../docs/AIEM_ROOT_CAUSE_ANALYSIS.md`

---

## 🔬 Technical Details

### Input Parameters
- **theta_i, theta_s** - Incident/scattered angles (degrees)
- **phi_s** - Scattered azimuth (degrees), use 180 for backscatter
- **kl** - Normalized correlation length (k × L)
- **ks** - Normalized RMS height (k × σ)
- **err, eri** - Real and imaginary parts of permittivity
- **itype** - Correlation type (1=Gaussian, 2=Exponential, 3=1.5-power)

### Output
- **VV, HH, HV, VH** - Backscatter coefficients (dB)

### Validity Range
- **ks:** 0.1 to 3.0 (optimal: 0.2 to 2.0)
- **kl:** 1.0 to 20.0 (optimal: 2.0 to 15.0)
- **theta:** 10° to 70° (optimal: 20° to 60°)

---

## 🔗 Related Files

### In This Directory
```
matlab/aiem_corrected/
├── *.m files (implementation)
└── *.md files (documentation)
```

### In Parent Project
```
docs/
├── AIEM_MATLAB_BUG_REPORT.md      ← Original bug report
├── AIEM_BUG_AUDIT_RESULTS.md      ← Bug audit
├── AIEM_ROOT_CAUSE_ANALYSIS.md    ← Root cause of remaining bias
└── AIEM_BUG_FIXES_SUMMARY.md      ← Executive summary

src/mwrtms/scattering/iem/
├── aiem.py                         ← Python reference implementation
├── fresnel_utils.py                ← Python Fresnel (with fixes)
├── transition.py                   ← Python transition (with fixes)
└── spectrum_aiem.py                ← Python spectrum (with fixes)
```

---

## ✅ Verification Checklist

### Installation
- [ ] Files copied to MATLAB directory
- [ ] Directory added to MATLAB path
- [ ] `test_aiem_corrected.m` runs without errors

### Functionality
- [ ] Basic computation works
- [ ] All correlation types work (1, 2, 3)
- [ ] Results are reasonable (-40 to 0 dB range)
- [ ] Monostatic reciprocity: HV ≈ VH

### Bug Fixes
- [ ] Fresnel: |R| ≤ 1 for lossy soils
- [ ] Normal incidence: rv0 = rh0
- [ ] Spectrum scaling laws verified
- [ ] No NaN or Inf in co-pol (VV, HH)

---

## 📞 Support & Resources

### Documentation
- **This directory:** Complete MATLAB implementation
- **`docs/` directory:** Analysis and bug reports
- **Python code:** Reference implementation

### Key Papers
1. Chen et al. (2003) - Original AIEM paper
2. Wu & Fung (1992) - Transition function
3. Fung et al. (1992) - IEM foundation

### Contact
- See main project README for contact information
- Check GitHub issues for known problems
- Review test scripts for usage examples

---

## 🎓 Learning Path

### Beginner
1. Read **README.md** - Understand what AIEM does
2. Run **test_aiem_corrected.m** - See it in action
3. Try **INSTALLATION.md** examples - Learn basic usage

### Intermediate
1. Study **BUG_FIXES.md** - Understand the corrections
2. Review **aiem_single_scattering.m** - See the algorithm
3. Run **compare_with_nmm3d.m** - Validate performance

### Advanced
1. Compare with Python implementation - Cross-reference
2. Read original papers - Theoretical foundation
3. Implement new transition function - Improve accuracy

---

## 📝 Version History

### Version 1.0 (2024)
- ✅ All 5 bugs fixed
- ✅ Complete documentation
- ✅ Test suite included
- ✅ NMM3D comparison script
- ⚠️ Legacy transition function (known +3-5 dB bias)

### Future (Planned)
- [ ] New S_p/S_p^(0) transition function
- [ ] Multiple scattering module
- [ ] Shadowing function
- [ ] Anisotropic correlation support

---

## 🏆 Status

**Implementation:** ✅ Complete  
**Bug Fixes:** ✅ All applied (5/5)  
**Testing:** ✅ Comprehensive  
**Documentation:** ✅ Complete  
**Validation:** ✅ Against NMM3D  

**Known Limitation:** ⚠️ +3-5 dB bias (legacy transition function)  
**Recommended Action:** Implement new transition method for <1 dB RMSE

---

## 🎯 Bottom Line

This is a **corrected, physically sound** implementation of AIEM with all identified bugs fixed. It's ready for use, with the caveat that the legacy transition function causes a systematic +3-5 dB bias vs NMM3D. For production applications requiring <1 dB RMSE, implement the new S_p/S_p^(0) transition method described in the bug report.

**Use this implementation if:**
- ✅ You need corrected AIEM (vs buggy original)
- ✅ You can tolerate +3-5 dB bias
- ✅ You need relative comparisons (bias cancels)
- ✅ You want a reference for Python implementation

**Upgrade to new transition if:**
- ⚠️ You need <1 dB RMSE vs NMM3D
- ⚠️ You need absolute calibration
- ⚠️ You're doing quantitative retrieval

---

**Last Updated:** 2024  
**Maintainer:** See main project  
**License:** Same as parent project
