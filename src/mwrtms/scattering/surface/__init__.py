"""Surface scattering models."""

from .base import SurfaceScattering
from .spm import SPMModel

# Legacy imports for backward compatibility
# AIEM and I2EM have been moved to scattering.iem
from ..iem.aiem import AIEMModel
from ..iem.i2em import I2EMModel

__all__ = ["SurfaceScattering", "SPMModel", "AIEMModel", "I2EMModel"]
