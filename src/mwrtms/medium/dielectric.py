"""Dielectric tensor utilities."""

from __future__ import annotations

import numpy as np

__all__ = ["DielectricTensor"]


class DielectricTensor:
    """Simple anisotropic dielectric tensor with encapsulated coefficients."""

    __slots__ = ("_eps_xx", "_eps_yy", "_eps_zz")
    _mutable_slots: set[str] = set()

    def __setattr__(self, name, value):
        if name in self.__slots__ and hasattr(self, name) and name not in self._mutable_slots:
            raise AttributeError(f"{self.__class__.__name__} attribute {name!r} is read-only")
        super().__setattr__(name, value)

    def __init__(self, eps_xx: complex, eps_yy: complex, eps_zz: complex) -> None:
        self._eps_xx = complex(eps_xx)
        self._eps_yy = complex(eps_yy)
        self._eps_zz = complex(eps_zz)

    @classmethod
    def isotropic(cls, eps: complex) -> "DielectricTensor":
        return cls(eps, eps, eps)

    @property
    def eps_xx(self) -> complex:
        return self._eps_xx

    @property
    def eps_yy(self) -> complex:
        return self._eps_yy

    @property
    def eps_zz(self) -> complex:
        return self._eps_zz

    @property
    def is_isotropic(self) -> bool:
        return bool(
            np.isclose(self._eps_xx, self._eps_yy)
            and np.isclose(self._eps_yy, self._eps_zz)
        )
