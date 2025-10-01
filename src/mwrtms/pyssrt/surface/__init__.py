from .aiem import AIEM
from .prism1 import PRISM1
from .smart import SMART
from .dubois95 import Dubois95
from .spm import SPM3D
from .i2em import *

__all__ = [
    "AIEM",
    "PRISM1",
    "SMART",
    "Dubois95",
    "SPM3D",
] + [name for name in globals() if name.startswith('I2EM')]
