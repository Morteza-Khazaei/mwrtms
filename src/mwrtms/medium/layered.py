"""Layered medium abstraction."""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

from .base import Medium

__all__ = ["Layer", "LayeredMedium"]


class Layer:
    """Single layer descriptor with encapsulated members."""

    __slots__ = ("_thickness", "_medium")
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, thickness: float | None, medium: Medium) -> None:
        if thickness is not None and thickness <= 0:
            raise ValueError("Layer thickness must be positive or None for semi-infinite")
        if not isinstance(medium, Medium):
            raise TypeError("medium must implement Medium")
        self._thickness = thickness
        self._medium = medium

    @property
    def thickness(self) -> float | None:
        return self._thickness

    @property
    def medium(self) -> Medium:
        return self._medium


class LayeredMedium:
    """Representation of a stratified medium."""

    __slots__ = ("_layers",)
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, layers: Iterable[Tuple[float | None, Medium]]) -> None:
        parsed: List[Layer] = []
        for thickness, medium in layers:
            parsed.append(Layer(thickness, medium))
        if not parsed:
            raise ValueError("At least one layer must be provided")
        self._layers: Sequence[Layer] = tuple(parsed)

    @property
    def layers(self) -> Sequence[Layer]:
        return self._layers

    def get_layer(self, depth: float) -> Layer:
        if depth < 0:
            raise ValueError("Depth must be non-negative")
        cumulative = 0.0
        for layer in self._layers:
            if layer.thickness is None:
                return layer
            cumulative += layer.thickness
            if depth < cumulative:
                return layer
        return self._layers[-1]

    def total_optical_depth(self, extinction_coefficients: Sequence[float]) -> float:
        if len(extinction_coefficients) != len(self._layers):
            raise ValueError("extinction_coefficients length must match number of layers")
        tau = 0.0
        for coeff, layer in zip(extinction_coefficients, self._layers):
            if layer.thickness is None:
                raise ValueError("Cannot compute optical depth with semi-infinite final layer")
            tau += coeff * layer.thickness
        return tau
