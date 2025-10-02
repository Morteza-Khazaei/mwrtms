"""Logging configuration helpers."""

from __future__ import annotations

import logging
from dataclasses import dataclass

__all__ = ["LoggingConfig"]


@dataclass
class LoggingConfig:
    """Simple logging configuration container."""

    level: int = logging.INFO
    console_output: bool = True

    @classmethod
    def development(cls) -> "LoggingConfig":
        return cls(level=logging.DEBUG, console_output=True)
