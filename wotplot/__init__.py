"A small library for creating and visualizing dot plot matrices."

from ._matrix import DotPlotMatrix
from ._make import MATCH, FWD, REV, BOTH
from ._viz import viz_spy, viz_imshow
from ._version import __version__

__all__ = [
    "DotPlotMatrix",
    "viz_spy",
    "viz_imshow",
    "style_viz_ax",
    "MATCH",
    "FWD",
    "REV",
    "BOTH",
    "__version__",
]
