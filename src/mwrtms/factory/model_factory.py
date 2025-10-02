"""Factory for creating scattering model instances."""

from __future__ import annotations

from typing import Dict, Type, TYPE_CHECKING

__all__ = ["ScatteringModelFactory"]

if TYPE_CHECKING:
    from ..core.radar_modes import RadarConfiguration
    from ..core.wave import ElectromagneticWave
    from ..scattering.base import ScatteringMechanism


class ScatteringModelFactory:
    """Factory maintaining a registry of scattering mechanisms."""

    _registry: Dict[str, Type] = {}

    @classmethod
    def register(cls, name: str, model_class: Type) -> None:
        cls._registry[name.lower()] = model_class

    @classmethod
    def create(cls, model_name: str, **kwargs):
        model_class = cls._registry.get(model_name.lower())
        if model_class is None:
            available = ", ".join(sorted(cls._registry))
            raise ValueError(f"Unknown model: {model_name}. Available: {available}")
        return model_class(**kwargs)

    @classmethod
    def list_models(cls) -> list[str]:
        return sorted(cls._registry)

    @classmethod
    def create_with_radar_config(
        cls,
        name: str,
        config: "RadarConfiguration",
        wave: "ElectromagneticWave",
        **kwargs,
    ) -> "ScatteringMechanism":
        model_class = cls._registry.get(name.lower())
        if model_class is None:
            available = ", ".join(sorted(cls._registry))
            raise ValueError(f"Unknown model: {name}. Available: {available}")

        existing_geometry = kwargs.get("geometry")
        if existing_geometry is not None and not config.matches_geometry(existing_geometry):
            raise ValueError("Provided geometry does not match radar configuration")

        kwargs.setdefault("geometry", config.geometry)
        kwargs.setdefault("wave", wave)
        return model_class(**kwargs)
