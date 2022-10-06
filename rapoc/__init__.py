from .__about__ import (
    __author__,
    __commit__,
    __copyright__,
    __email__,
    __license__,
    __summary__,
    __title__,
    __pkg_name__,
    __url__,
    __branch__,
    __base_dir__
)
from .__version__ import __version__

__all__ = [
    "__author__",
    "__commit__",
    "__copyright__",
    "__email__",
    "__license__",
    "__summary__",
    "__title__",
    "__pkg_name__",
    "__url__",
    "__version__",
    "__branch__",
    "__base_dir__"
]

from . import models
from .models import Rosseland, Planck, Rayleigh
