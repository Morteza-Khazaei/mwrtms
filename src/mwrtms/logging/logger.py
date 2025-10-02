"""Centralised logger management."""

from __future__ import annotations

import logging

__all__ = ["MWLogger"]


class MWLogger:
    """Singleton-style access to configured loggers."""

    _loggers: dict[str, logging.Logger] = {}

    @classmethod
    def get_logger(cls, name: str) -> logging.Logger:
        if name not in cls._loggers:
            logger = logging.getLogger(name)
            logger.setLevel(logging.INFO)
            if not logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
                handler.setFormatter(formatter)
                logger.addHandler(handler)
            cls._loggers[name] = logger
        return cls._loggers[name]
