import cv2
import numpy as np
import matplotlib
from math import ceil
from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator

# Based on the opencv tutorial linked in viz_binary()'s docstring
DILATION_KERNEL = np.ones((5, 5), np.uint8)


def _tidy_dotplot_viz_ax(ax, m):
    # Force matplotlib to only use integer ticks: https://stackoverflow.com/a/38096332
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Hide the ticks / tick labels. 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Based on the mutation matrix drawn in
    # https://nbviewer.org/github/fedarko/strainFlye/blob/main/docs/SheepGutExample.ipynb
    ax.tick_params(
        top=True, bottom=True, left=True, right=True,
        labeltop=True, labelbottom=True, labelleft=True, labelright=True
    )

    ax.set_xlabel(f"$s_1$ ({len(m.s1):,} nt)", fontsize=18)
    ax.set_ylabel(f"$s_2$ ({len(m.s2):,} nt)", fontsize=18)

def viz_binary(m, num_dilation_iterations="auto", title=None, ax=None):
    """Visualizes a DotPlotMatrix object.

    This is intended to be used on binary DotPlotMatrix objects only.

    References
    ----------
    https://docs.opencv.org/4.x/d9/d61/tutorial_py_morphological_ops.html
        OpenCV tutorial on dilation and other morphological operations.

    https://stackoverflow.com/q/44618675
        The answers to this question showed me the light (of image processing
        101) (I didn't know what a dilation was)
    """
    if not m.binary:
        raise ValueError(
            "Can't use viz_binary() on a DotPlotMatrix that is not binary."
        )

    if num_dilation_iterations == "auto":
        maxlen = max(len(m.s1), len(m.s2))
        # Don't bother doing dilation for small matrices
        if maxlen < 500:
            num_dilation_iterations = 0
        else:
            # sloppy stuff, todo replace (also check stuff in [500, 2500))
            num_dilation_iterations = 4 + ceil(maxlen - 2500) / 2500

    if num_dilation_iterations > 0:
        matrix_to_show = cv2.dilate(m.mat, DILATION_KERNEL, iterations=num_dilation_iterations)
    else:
        matrix_to_show = m.mat

    return_mplobjs = False
    if ax is None:
        return_mplobjs = True
        fig, ax = pyplot.subplots()

    ax.imshow(matrix_to_show, cmap="gray_r")

    _tidy_dotplot_viz_ax(ax, m)

    if title is not None:
        ax.set_title(title, fontsize=18)

    if return_mplobjs:
        return fig, ax
