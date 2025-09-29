"""Build helpers for pySSRT extension modules."""

from __future__ import annotations

from setuptools.command.build_ext import build_ext


class BuildExt(build_ext):
    """Custom build_ext that cythonizes sources and adds NumPy headers."""

    _already_cythonized = False

    def finalize_options(self):  # noqa: D401 - inherited docstring from distutils
        super().finalize_options()
        self._ensure_numpy_include()
        self._cythonize_extensions()

    def _ensure_numpy_include(self) -> None:
        import numpy  # noqa: WPS433 (delayed import during build)

        include_dir = numpy.get_include()
        for extension in self.extensions:
            if include_dir not in extension.include_dirs:
                extension.include_dirs.append(include_dir)

    def _cythonize_extensions(self) -> None:
        if self._already_cythonized:
            return
        try:
            from Cython.Build import cythonize
        except ImportError as exc:  # pragma: no cover - build-time environment
            raise RuntimeError(
                "Cython is required to build the optional _i2em_cy extension. "
                "Install with the 'build' extras or ensure Cython is available."
            ) from exc

        self.extensions = cythonize(  # type: ignore[assignment]
            self.extensions,
            compiler_directives={"language_level": "3"},
        )
        self._already_cythonized = True


__all__ = ["BuildExt"]
