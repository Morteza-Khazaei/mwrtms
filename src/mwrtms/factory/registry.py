"""Registration decorator for scattering models."""

from __future__ import annotations

from typing import TypeVar

T = TypeVar("T")


def register_model(name: str):
    """Decorator that auto-registers scattering models with the factory."""

    def decorator(cls: T) -> T:
        from .model_factory import ScatteringModelFactory

        ScatteringModelFactory.register(name, cls)
        return cls

    return decorator
