import numpy as np

from mwrtms.medium.surface import (
    Surface,
    SurfaceAnalyzer,
    SyntheticSurfaceGenerator,
)


def test_surface_creation() -> None:
    rng = np.random.default_rng(0)
    heights = rng.normal(scale=1e-6, size=(100, 120))
    surface = Surface(heights, (1e-3, 1.2e-3))

    assert surface.shape == (100, 120)
    assert surface.physical_size == (1e-3, 1.2e-3)


def test_rms_height_positive() -> None:
    rng = np.random.default_rng(1)
    heights = rng.normal(scale=1e-6, size=(64, 64))
    surface = Surface(heights, (1e-3, 1e-3))
    assert surface.rms_height() > 0.0


def test_correlation_length_reasonable() -> None:
    generator = SyntheticSurfaceGenerator((128, 128), (1e-3, 1e-3)).setRandomSeed(123)
    surface = generator.generate()
    corr_length = surface.correlation_length()
    assert 0.0 < corr_length < 1e-3


def test_isotropy_detection() -> None:
    generator = SyntheticSurfaceGenerator((128, 128), (1e-3, 1e-3)).setRandomSeed(42)
    surface = generator.generate()
    assert surface.is_isotropic()


def test_surface_analyzer_outputs_keys() -> None:
    generator = SyntheticSurfaceGenerator((64, 64), (5e-4, 5e-4)).setRandomSeed(10)
    surface = generator.generate()
    analysis = SurfaceAnalyzer.analyze(surface)
    expected_keys = {
        "rms_height",
        "mean_height",
        "peak_to_valley",
        "rms_slope",
        "is_isotropic",
        "anisotropy_ratio",
        "correlation_length",
        "correlation_length_spectral",
    }
    assert expected_keys.issubset(analysis.keys())


def test_numpy_integration_helpers() -> None:
    generator = SyntheticSurfaceGenerator((32, 32), (1e-3, 1e-3)).setRandomSeed(99)
    surface = generator.generate()

    heights = np.array(surface)
    assert heights.shape == surface.shape

    subset = surface[0:5, 0:5]
    assert subset.shape == (5, 5)

    mean_val = np.mean(surface)
    assert isinstance(mean_val, (float, np.floating))


def test_with_physical_size_preserves_heights() -> None:
    generator = SyntheticSurfaceGenerator((32, 32), (1e-3, 1e-3)).setRandomSeed(5)
    surface = generator.generate()
    resized = surface.with_physical_size((2e-3, 1.5e-3))

    assert resized.shape == surface.shape
    np.testing.assert_allclose(resized.heights, surface.heights)
    assert resized.physical_size == (2e-3, 1.5e-3)


def test_analytic_spectrum_positive() -> None:
    generator = SyntheticSurfaceGenerator((64, 64), (1e-3, 1e-3)).setRandomSeed(21)
    surface = generator.generate()
    value = surface.analytic_spectrum(0, 0.0, 0.0)
    assert value >= 0.0


def test_builder_matches_targets() -> None:
    target_rms = 2e-3
    target_corr = 8e-3
    surface = build_surface_from_statistics(target_rms, target_corr, correlation_type="gaussian")

    assert abs(surface.rms_height() - target_rms) / target_rms < 1e-6
    assert surface.metadata.get("correlation_type") == "gaussian"
    assert abs(surface.correlation_length() - target_corr) < 1e-9
