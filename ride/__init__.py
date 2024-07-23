from .ride import ride_call
from .cfg import RideCfg
from .correct import correct_trials

# Get current package version
try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")
