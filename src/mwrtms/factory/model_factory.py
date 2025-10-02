"""Factory for creating scattering model instances."""

from __future__ import annotations

from typing import Dict, Type

__all__ = ["ScatteringModelFactory"]


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
