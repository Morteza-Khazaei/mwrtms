"""Factory utilities for scattering models."""

from __future__ import annotations

from typing import Dict, Type

from ..scattering import (
    AIEMModel,
    Dubois95Model,
    I2EMModel,
    IEMModel,
    PRISM1Model,
    SMARTModel,
    SPM3DModel,
    SSRTModel,
    ScatteringMechanism,
)

__all__ = ["ScatteringModelFactory", "register_model"]


class ScatteringModelFactory:
    """Factory maintaining a registry of scattering mechanisms."""

    _registry: Dict[str, Type[ScatteringMechanism]] = {
        "aiem": AIEMModel,
        "i2em": I2EMModel,
        "iem": IEMModel,
        "prism1": PRISM1Model,
        "spm3d": SPM3DModel,
        "smart": SMARTModel,
        "dubois95": Dubois95Model,
        "ssrt": SSRTModel,
    }

    @classmethod
    def register(cls, name: str, model_class: Type[ScatteringMechanism]) -> None:
        cls._registry[name.lower()] = model_class

    @classmethod
    def create(cls, name: str, **kwargs) -> ScatteringMechanism:
        key = name.lower()
        if key not in cls._registry:
            available = ", ".join(sorted(cls._registry))
            raise ValueError(f"Unknown model '{name}'. Available models: {available}")
        return cls._registry[key](**kwargs)

    @classmethod
    def list_models(cls) -> list[str]:
        return sorted(cls._registry)


def register_model(name: str):
    """Decorator for registering additional models."""

    def decorator(model_class: Type[ScatteringMechanism]) -> Type[ScatteringMechanism]:
        ScatteringModelFactory.register(name, model_class)
        return model_class

    return decorator
