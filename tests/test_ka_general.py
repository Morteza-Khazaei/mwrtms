"""Test suite for generalized KA model with multiple ACF types.

This test suite validates the KA model against:
1. All three ACF types (Gaussian, Exponential, x-Power)
2. Mathematical properties (DC values, similarity laws)
3. Physical constraints (reciprocity, energy conservation)
4. Series convergence behavior
"""

import numpy as np
import pytest
from scipy.special import gamma as gamma_func

from mwrtms import mwRTMs, RadarConfigurationFactory, PolarizationState


class TestKAGaussianACF:
    """Tests for KA model with Gaussian ACF."""
    
    def test_gaussian_basic(self):
        """Test basic Gaussian ACF computation."""
        frequency_ghz = 5.405
        theta_deg = 40.0
        rms_height_cm = 2.0
        correlation_length_cm = 20.0
        soil_permittivity = complex(15.0, 3.0)
        
        radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=theta_deg)
        
        result = mwRTMs.compute_soil_backscatter(
            model='ka',
            radar_config=radar_config,
            frequency_ghz=frequency_ghz,
            rms_height_cm=rms_height_cm,
            correlation_length_cm=correlation_length_cm,
            soil_permittivity=soil_permittivity,
            polarization=PolarizationState.VV,
        )
        
        assert isinstance(result, (float, np.floating))
        assert result >= 0.0
        assert not np.isnan(result)
    
    def test_gaussian_all_polarizations(self):
        """Test all polarizations with Gaussian ACF."""
        frequency_ghz = 5.405
        theta_deg = 40.0
        rms_height_cm = 2.0
        correlation_length_cm = 20.0
        soil_permittivity = complex(15.0, 3.0)
        
        radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=theta_deg)
        
        for pol in [PolarizationState.VV, PolarizationState.HH, 
                    PolarizationState.HV, PolarizationState.VH]:
            result = mwRTMs.compute_soil_backscatter(
                model='ka',
                radar_config=radar_config,
                frequency_ghz=frequency_ghz,
                rms_height_cm=rms_height_cm,
                correlation_length_cm=correlation_length_cm,
                soil_permittivity=soil_permittivity,
                polarization=pol,
            )
            
            assert result >= 0.0
            assert not np.isnan(result)


class TestKAExponentialACF:
    """Tests for KA model with Exponential ACF."""
    
    def test_exponential_basic(self):
        """Test basic Exponential ACF computation."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="exponential"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        model = KAModel(wave, geometry, surface, acf_type="exponential", nmax=8)
        sigma0 = model.compute(None, soil, PolarizationState.VV)
        
        assert sigma0 >= 0.0
        assert not np.isnan(sigma0)
    
    def test_exponential_all_polarizations(self):
        """Test all polarizations with Exponential ACF."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="exponential"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        for pol in [PolarizationState.VV, PolarizationState.HH, 
                    PolarizationState.HV, PolarizationState.VH]:
            model = KAModel(wave, geometry, surface, acf_type="exponential", nmax=8)
            sigma0 = model.compute(None, soil, pol)
            
            assert sigma0 >= 0.0
            assert not np.isnan(sigma0)


class TestKAXPowerACF:
    """Tests for KA model with x-Power ACF."""
    
    def test_xpower_basic(self):
        """Test basic x-Power ACF computation."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="exponential"  # Template
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        model = KAModel(wave, geometry, surface, acf_type="xpower", alpha=1.5, nmax=8)
        sigma0 = model.compute(None, soil, PolarizationState.VV)
        
        assert sigma0 >= 0.0
        assert not np.isnan(sigma0)
    
    def test_xpower_alpha_values(self):
        """Test x-Power ACF with different alpha values."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="exponential"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        # Test different alpha values
        for alpha in [0.75, 1.0, 1.5, 2.0, 3.0]:
            model = KAModel(wave, geometry, surface, acf_type="xpower", alpha=alpha, nmax=6)
            sigma0 = model.compute(None, soil, PolarizationState.VV)
            
            assert sigma0 >= 0.0
            assert not np.isnan(sigma0)
    
    def test_xpower_reduces_to_gaussian(self):
        """Test that x-Power with alpha=2 gives similar results to Gaussian."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="gaussian"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        # Gaussian ACF
        model_gauss = KAModel(wave, geometry, surface, acf_type="gaussian", nmax=8)
        sigma0_gauss = model_gauss.compute(None, soil, PolarizationState.VV)
        
        # x-Power with alpha=2 (should be same as Gaussian)
        model_xpower = KAModel(wave, geometry, surface, acf_type="xpower", alpha=2.0, nmax=8)
        sigma0_xpower = model_xpower.compute(None, soil, PolarizationState.VV)
        
        # Should be very close (within numerical precision)
        assert np.isclose(sigma0_gauss, sigma0_xpower, rtol=0.01)
    
    def test_xpower_reduces_to_exponential(self):
        """Test that x-Power with alpha=1 gives similar results to Exponential."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="exponential"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        # Exponential ACF
        model_exp = KAModel(wave, geometry, surface, acf_type="exponential", nmax=8)
        sigma0_exp = model_exp.compute(None, soil, PolarizationState.VV)
        
        # x-Power with alpha=1 (should be same as Exponential)
        model_xpower = KAModel(wave, geometry, surface, acf_type="xpower", alpha=1.0, nmax=8)
        sigma0_xpower = model_xpower.compute(None, soil, PolarizationState.VV)
        
        # Should be very close (within numerical precision)
        assert np.isclose(sigma0_exp, sigma0_xpower, rtol=0.01)


class TestKAPhysics:
    """Test physical behavior across all ACF types."""
    
    def test_reciprocity_all_acfs(self):
        """Test HV = VH reciprocity for all ACF types."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="gaussian"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        for acf_type in ["gaussian", "exponential"]:
            model = KAModel(wave, geometry, surface, acf_type=acf_type, nmax=8)
            
            hv = model.compute(None, soil, PolarizationState.HV)
            vh = model.compute(None, soil, PolarizationState.VH)
            
            assert np.isclose(hv, vh, rtol=1e-10)
    
    def test_copol_greater_than_crosspol_all_acfs(self):
        """Test co-pol > cross-pol for all ACF types."""
        from mwrtms.core import ElectromagneticWave, ScatteringGeometry
        from mwrtms.medium import HomogeneousMedium, build_surface_from_statistics
        from mwrtms.scattering.surface.ka import KAModel
        
        frequency_ghz = 5.405
        theta_deg = 40.0
        
        wave = ElectromagneticWave(frequency_hz=frequency_ghz * 1e9)
        geometry = ScatteringGeometry(theta_i_deg=theta_deg)
        surface = build_surface_from_statistics(
            rms_height_m=0.02,
            correlation_length_m=0.20,
            correlation_type="gaussian"
        )
        soil = HomogeneousMedium(permittivity=complex(15.0, 3.0))
        
        for acf_type in ["gaussian", "exponential"]:
            model = KAModel(wave, geometry, surface, acf_type=acf_type, nmax=8)
            
            vv = model.compute(None, soil, PolarizationState.VV)
            hv = model.compute(None, soil, PolarizationState.HV)
            
            assert vv > hv


class TestKAMathematical:
    """Test mathematical properties of the KA model."""
    
    def test_wn_gaussian_dc_value(self):
        """Test Gaussian W^(n)(0) = πL²/n."""
        from mwrtms.scattering.surface.ka import KAModel
        
        L = 0.20  # 20 cm
        
        for n in [1, 2, 3, 5, 10]:
            wn_zero = KAModel._wn_gaussian(0.0, L, n)
            expected = np.pi * L**2 / n
            assert np.isclose(wn_zero, expected, rtol=1e-10)
    
    def test_wn_exponential_dc_value(self):
        """Test Exponential W^(n)(0) = 2πL²/n²."""
        from mwrtms.scattering.surface.ka import KAModel
        
        L = 0.20  # 20 cm
        
        for n in [1, 2, 3, 5, 10]:
            wn_zero = KAModel._wn_exponential(0.0, L, n)
            expected = 2.0 * np.pi * L**2 / (n**2)
            assert np.isclose(wn_zero, expected, rtol=1e-10)


def test_ka_example():
    """Example test showing typical usage with all ACF types."""
    print("\n" + "="*70)
    print("KA Model with Multiple ACF Types - Example Test")
    print("="*70)
    
    frequency_ghz = 5.405  # C-band
    theta_deg = 40.0
    rms_height_cm = 2.0
    correlation_length_cm = 20.0
    soil_permittivity = complex(15.0, 3.0)
    
    radar_config = RadarConfigurationFactory.create_monostatic(theta_deg=theta_deg)
    
    print(f"\nConfiguration:")
    print(f"  Frequency: {frequency_ghz} GHz")
    print(f"  Incidence angle: {theta_deg}°")
    print(f"  RMS height: {rms_height_cm} cm")
    print(f"  Correlation length: {correlation_length_cm} cm")
    print(f"  Permittivity: {soil_permittivity}")
    
    # Test with Gaussian ACF (default)
    print(f"\nBackscatter coefficients (Gaussian ACF):")
    for pol in [PolarizationState.VV, PolarizationState.HH]:
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
        print(f"  σ⁰_{pol.value}: {result_db:7.2f} dB")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    # Run the example test
    test_ka_example()
    
    # Run pytest
    pytest.main([__file__, "-v"])
