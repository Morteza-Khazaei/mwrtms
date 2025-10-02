"""Facade providing simplified APIs."""

from .api import mwRTMs

mwRTMsFacade = mwRTMs

__all__ = ["mwRTMs", "mwRTMsFacade"]
