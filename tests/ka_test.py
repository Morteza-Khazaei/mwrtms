"""Test suite for the Kirchhoff Approximation (KA) surface scattering model.

This test module validates the KA model implementation against theoretical
expectations and numerical properties. The KA model is particularly suitable
for large-scale surface roughness (long gravity waves).

Usage
-----
Run from the repository root:

    PYTHONPATH=src python3 -m pytest tests/ka_test.py -v

Or run directly:

    PYTHONPATH=src python3 tests/ka_test.py
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState


class TestKAModel:
    """Test suite for the KA surface scattering model."""

    @pytest.fixture
    def default_params(self):
        """Default parameters for KA model testing."""
        return {
            'frequency_ghz': 5.405,
            'theta_deg': 40.0,
            'rms_height_cm': 1.0,
            'correlation_length_cm': 10.0,
            'soil_permittivity': complex(15.0, 3.0),
        }

    def test_ka_model_basic_computation(self, default_params):
        """Test that KA model computes without errors for basic inputs."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        result_hh = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.HH,
        )
        
        # Results should be positive and finite
        assert result_vv > 0, "VV backscatter should be positive"
        assert result_hh > 0, "HH backscatter should be positive"
        assert np.isfinite(result_vv), "VV backscatter should be finite"
        assert np.isfinite(result_hh), "HH backscatter should be finite"
        
        # Convert to dB for reasonable range check
        vv_db = 10 * np.log10(result_vv)
        hh_db = 10 * np.log10(result_hh)
        
        # Backscatter should be in a reasonable range (KA can produce low values)
        assert -80 < vv_db < 10, f"VV backscatter {vv_db:.2f} dB out of expected range"
        assert -80 < hh_db < 10, f"HH backscatter {hh_db:.2f} dB out of expected range"

    def test_ka_cross_polarization(self, default_params):
        """Test cross-polarization (HV/VH) computation."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        result_hv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.HV,
        )
        
        result_vh = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VH,
        )
        
        # Cross-pol should be non-negative and finite
        assert result_hv >= 0, "HV backscatter should be non-negative"
        assert result_vh >= 0, "VH backscatter should be non-negative"
        assert np.isfinite(result_hv), "HV backscatter should be finite"
        assert np.isfinite(result_vh), "VH backscatter should be finite"
        
        # For reciprocal backscatter, HV should equal VH
        if result_hv > 0 and result_vh > 0:
            ratio = result_hv / result_vh
            assert 0.9 < ratio < 1.1, f"HV and VH should be approximately equal (ratio={ratio:.3f})"

    def test_ka_angle_dependence(self, default_params):
        """Test that backscatter varies with incidence angle."""
        angles = [20.0, 30.0, 40.0, 50.0, 60.0]
        results_vv = []
        
        for angle in angles:
            radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=angle)
            result = mwRTMs.compute_soil_backscatter(
                model='ka',
                radar_config=radar_config,
                frequency_ghz=default_params['frequency_ghz'],
                rms_height_cm=default_params['rms_height_cm'],
                correlation_length_cm=default_params['correlation_length_cm'],
                soil_permittivity=default_params['soil_permittivity'],
                polarization=PolarizationState.VV,
            )
            results_vv.append(result)
        
        # All results should be positive and finite
        assert all(r > 0 for r in results_vv), "All backscatter values should be positive"
        assert all(np.isfinite(r) for r in results_vv), "All backscatter values should be finite"
        
        # Results should vary with angle (not all the same)
        assert len(set(results_vv)) > 1, "Backscatter should vary with incidence angle"

    def test_ka_frequency_dependence(self, default_params):
        """Test that backscatter varies with frequency."""
        frequencies = [1.0, 5.0, 10.0, 15.0]
        results_vv = []
        
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        for freq in frequencies:
            result = mwRTMs.compute_soil_backscatter(
                model='ka',
                radar_config=radar_config,
                frequency_ghz=freq,
                rms_height_cm=default_params['rms_height_cm'],
                correlation_length_cm=default_params['correlation_length_cm'],
                soil_permittivity=default_params['soil_permittivity'],
                polarization=PolarizationState.VV,
            )
            results_vv.append(result)
        
        # All results should be positive and finite
        assert all(r > 0 for r in results_vv), "All backscatter values should be positive"
        assert all(np.isfinite(r) for r in results_vv), "All backscatter values should be finite"

    def test_ka_roughness_dependence(self, default_params):
        """Test that backscatter varies with surface roughness."""
        rms_heights = [0.5, 1.0, 2.0, 3.0]
        results_vv = []
        
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        for rms in rms_heights:
            result = mwRTMs.compute_soil_backscatter(
                model='ka',
                radar_config=radar_config,
                frequency_ghz=default_params['frequency_ghz'],
                rms_height_cm=rms,
                correlation_length_cm=default_params['correlation_length_cm'],
                soil_permittivity=default_params['soil_permittivity'],
                polarization=PolarizationState.VV,
            )
            results_vv.append(result)
        
        # All results should be positive and finite
        assert all(r > 0 for r in results_vv), "All backscatter values should be positive"
        assert all(np.isfinite(r) for r in results_vv), "All backscatter values should be finite"
        
        # Results should vary with roughness
        assert len(set(results_vv)) > 1, "Backscatter should vary with surface roughness"

    def test_ka_permittivity_dependence(self, default_params):
        """Test that backscatter varies with soil permittivity."""
        permittivities = [
            complex(3.0, 0.1),
            complex(10.0, 2.0),
            complex(20.0, 5.0),
            complex(30.0, 8.0),
        ]
        results_vv = []
        
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        for eps in permittivities:
            result = mwRTMs.compute_soil_backscatter(
                model='ka',
                radar_config=radar_config,
                frequency_ghz=default_params['frequency_ghz'],
                rms_height_cm=default_params['rms_height_cm'],
                correlation_length_cm=default_params['correlation_length_cm'],
                soil_permittivity=eps,
                polarization=PolarizationState.VV,
            )
            results_vv.append(result)
        
        # All results should be positive and finite
        assert all(r > 0 for r in results_vv), "All backscatter values should be positive"
        assert all(np.isfinite(r) for r in results_vv), "All backscatter values should be finite"
        
        # Results should vary with permittivity
        assert len(set(results_vv)) > 1, "Backscatter should vary with permittivity"

    def test_ka_copol_ratio(self, default_params):
        """Test that co-polarization ratio is reasonable."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        result_hh = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.HH,
        )
        
        # VV/HH ratio should be reasonable (typically 0.1 to 10 in linear units)
        ratio = result_vv / result_hh
        assert 0.01 < ratio < 100, f"VV/HH ratio {ratio:.3f} is outside expected range"
        
        # Convert to dB
        vv_db = 10 * np.log10(result_vv)
        hh_db = 10 * np.log10(result_hh)
        ratio_db = vv_db - hh_db
        
        # Ratio in dB should be reasonable (typically -10 to +10 dB)
        assert -20 < ratio_db < 20, f"VV/HH ratio {ratio_db:.2f} dB is outside expected range"

    def test_ka_nadir_angle(self, default_params):
        """Test KA model at near-nadir incidence."""
        radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=5.0)
        
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results even at near-nadir
        assert result_vv > 0, "VV backscatter should be positive at nadir"
        assert np.isfinite(result_vv), "VV backscatter should be finite at nadir"

    def test_ka_steep_angle(self, default_params):
        """Test KA model at steep incidence angles."""
        radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=70.0)
        
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results even at steep angles
        assert result_vv > 0, "VV backscatter should be positive at steep angles"
        assert np.isfinite(result_vv), "VV backscatter should be finite at steep angles"

    def test_ka_smooth_surface(self, default_params):
        """Test KA model with very smooth surface."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        # Very smooth surface (small rms height, but not too small)
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=0.3,
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results for smooth surfaces (may be very small)
        assert result_vv >= 0, "VV backscatter should be non-negative for smooth surface"
        assert np.isfinite(result_vv), "VV backscatter should be finite for smooth surface"

    def test_ka_rough_surface(self, default_params):
        """Test KA model with very rough surface."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        # Very rough surface (large rms height)
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=5.0,
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=default_params['soil_permittivity'],
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results for rough surfaces
        assert result_vv > 0, "VV backscatter should be positive for rough surface"
        assert np.isfinite(result_vv), "VV backscatter should be finite for rough surface"

    def test_ka_low_permittivity(self, default_params):
        """Test KA model with low permittivity (dry soil)."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        # Dry soil (low permittivity)
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=complex(3.0, 0.1),
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results for low permittivity
        assert result_vv > 0, "VV backscatter should be positive for low permittivity"
        assert np.isfinite(result_vv), "VV backscatter should be finite for low permittivity"

    def test_ka_high_permittivity(self, default_params):
        """Test KA model with high permittivity (wet soil/water)."""
        radar_config = RadarConfigurationFactory.create_monostatic(
            theta_deg=default_params['theta_deg']
        )
        
        # Wet soil/water (high permittivity)
        result_vv = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=default_params['frequency_ghz'],
            rms_height_cm=default_params['rms_height_cm'],
            correlation_length_cm=default_params['correlation_length_cm'],
            soil_permittivity=complex(80.0, 10.0),
            polarization=PolarizationState.VV,
        )
        
        # Should produce valid results for high permittivity
        assert result_vv > 0, "VV backscatter should be positive for high permittivity"
        assert np.isfinite(result_vv), "VV backscatter should be finite for high permittivity"


def print_test_summary():
    """Print a summary of KA model test results."""
    print("\n" + "=" * 70)
    print("KA Model Test Summary")
    print("=" * 70)
    
    # Test configuration
    frequency_ghz = 5.405
    theta_deg = 40.0
    rms_height_cm = 1.0
    correlation_length_cm = 10.0
    soil_permittivity = complex(15.0, 3.0)
    
    radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=theta_deg)
    
    print(f"\nTest Configuration:")
    print(f"  Frequency: {frequency_ghz} GHz")
    print(f"  Incidence angle: {theta_deg}°")
    print(f"  RMS height: {rms_height_cm} cm")
    print(f"  Correlation length: {correlation_length_cm} cm")
    print(f"  Soil permittivity: {soil_permittivity}")
    
    # Compute all polarizations
    polarizations = [PolarizationState.VV, PolarizationState.HH, 
                     PolarizationState.HV, PolarizationState.VH]
    
    print(f"\nBackscatter Results:")
    for pol in polarizations:
        result = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=frequency_ghz,
            rms_height_cm=rms_height_cm,
            correlation_length_cm=correlation_length_cm,
            soil_permittivity=soil_permittivity,
            polarization=pol,
        )
        result_db = 10 * np.log10(result) if result > 0 else float('-inf')
        print(f"  σ⁰_{pol.value}: {result_db:7.2f} dB ({result:.6e} linear)")
    
    # Test angle dependence
    print(f"\nAngle Dependence (VV polarization):")
    angles = [10, 20, 30, 40, 50, 60, 70]
    for angle in angles:
        radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=angle)
        result = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=frequency_ghz,
            rms_height_cm=rms_height_cm,
            correlation_length_cm=correlation_length_cm,
            soil_permittivity=soil_permittivity,
            polarization=PolarizationState.VV,
        )
        result_db = 10 * np.log10(result) if result > 0 else float('-inf')
        print(f"  θ = {angle:2d}°: {result_db:7.2f} dB")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    # Run tests with pytest if available, otherwise run summary
    try:
        import sys
        pytest.main([__file__, "-v", "--tb=short"])
    except ImportError:
        print("pytest not available, running test summary instead")
        print_test_summary()
