__version__ = "0.0.0"  # Replaced by poetry-dynamic-versioning at build time

import pandas as _pd

# SADIE stores strings as object-dtype columns and treats missing values as float
# NaN (the AirrTable/NumberingResults logic relies on ``isinstance(x, float)``).
# pandas 3.0 (PDEP-14) defaults to a new StringDtype; opt back out so behavior
# matches pandas 2.x. The option does not exist before pandas 2.1, hence the guard.
try:
    _pd.set_option("future.infer_string", False)
except (KeyError, _pd.errors.OptionError):
    pass

from . import airr, numbering, receptor, reference, renumbering

__all__ = ["renumbering", "receptor", "airr", "reference", "numbering"]
