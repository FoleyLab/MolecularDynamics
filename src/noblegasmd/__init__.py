"""noblegasmd: a numba port of the NVE Lennard-Jones noble-gas MD engine
originally implemented in MD.cpp by Foley, Sweet, and Akinfenwa (GPLv3).

See https://doi.org/10.1021/acs.jchemed.7b00790 for the original article.
"""
from .api import RunResult, run, sweep
from .units import GAS_CONSTANTS, KB_SI, NA, GasConstants, get_gas_constants

__version__ = "0.1.0"

__all__ = [
    "run",
    "sweep",
    "RunResult",
    "GAS_CONSTANTS",
    "GasConstants",
    "get_gas_constants",
    "KB_SI",
    "NA",
    "__version__",
]
