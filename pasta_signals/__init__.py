"""
Old name of the pastatemp package, kept so existing scripts still run.
Use `import pastatemp` in new code.
"""

import sys
import warnings

import pastatemp as _pt
from pastatemp import *  # noqa: F401,F403
from pastatemp import (annual, core, coverage, diel, export, pairing, phase,  # noqa: F401
                       readers)

__version__ = _pt.__version__

for _m in ("annual", "core", "coverage", "diel", "export", "pairing", "phase", "readers"):
    sys.modules[__name__ + "." + _m] = sys.modules["pastatemp." + _m]

warnings.warn("pasta_signals is now called pastatemp; please `import pastatemp`.",
              DeprecationWarning, stacklevel=2)
