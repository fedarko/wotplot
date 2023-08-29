import scipy
from scipy.ndimage import binary_dilation
import numpy as np
from math import ceil
from matplotlib import pyplot

# Based on the opencv tutorial linked in viz_binary()'s docstring
DILATION_KERNEL = np.ones((5, 5), np.uint8)


def _draw_sparse_binary_matrix(ax, mat):
    # We don't need to worry about yorder here, since that should already have
    # been "applied" to the matrix
    ax.set_xlim(-0.5, mat.shape[1] - 0.5)
    ax.set_ylim(mat.shape[0] - 0.5, -0.5)

    ax.set_aspect("equal")
    # if we don't do this then the bbox is wrong
    ax.apply_aspect()

    # Figure out how big each point in the visualization should be.
    #
    # ax.get_window_extent() gives us the dimensions of *just* the interior of
    # the plot (ignoring ticks, etc.: https://stackoverflow.com/a/64355139),
    # and we can use this to figure out the edge lengths of each cell's
    # square in the plot.
    bbox = ax.get_window_extent()
    # We could also compute this using bbox.height / mat.shape[0] -- I've
    # tested it, and these are equivalent.
    cell_side_length = bbox.width / mat.shape[1]
    # The "s" parameter of ax.scatter() is in terms of point area. This is why
    # we need to square the side length in order to get a value of s that
    # *just* fills up each cell's "square".
    cell_area = cell_side_length**2

    # Insane sidenote: to illustrate the "cell squares" that I'm talking about
    # here, you can run the following lines after creating a visualization.
    # (Run them after _tidy_dotplot_viz_ax(), since otherwise that'll override
    # these ticks.)
    #
    # ax.set_xticks([x - 0.5 for x in range(0, mat.shape[1])])
    # ax.set_yticks([x - 0.5 for x in range(0, mat.shape[0])])
    # ax.grid()
    #
    # Note that the grid might look slightly misaligned from the points
    # depending on your resolution -- I was really confused about this for a
    # while, and then I tried adding the same ticks & grid lines for the output
    # of ax.imshow() on the same matrix but dense, and whaddaya know the grid
    # lines there were also slightly off (but in a different way). This was
    # also the case with betterspy (https://github.com/nschloe/betterspy).
    # In any case, I figure this current code is Good Enough (TM).

    # The .col and .row attributes of mat refer to the nonzero entries' cells.
    ax.scatter(mat.col, mat.row, marker="s", color="#000", s=cell_area)


def _tidy_dotplot_viz_ax(ax, m):
    # Hide the ticks / tick labels.
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Based on the mutation matrix drawn in
    # https://nbviewer.org/github/fedarko/strainFlye/blob/main/docs/SheepGutExample.ipynb
    ax.tick_params(
        top=True,
        bottom=True,
        left=True,
        right=True,
        labeltop=True,
        labelbottom=True,
        labelleft=True,
        labelright=True,
    )

    ax.set_xlabel(f"$s_1$ ({len(m.s1):,} nt) \u2192", fontsize=18)
    if m.yorder == "BT":
        # -->, which gets turned into an up arrow when we rotate the yax label
        yarr = "\u2192"
    elif m.yorder == "TB":
        # <--, which gets turned into a down arrow when we rotate the yax label
        yarr = "\u2190"
    else:
        raise ValueError(f"Unrecognized yorder?: {m.yorder}")
    ax.set_ylabel(f"$s_2$ ({len(m.s2):,} nt) {yarr}", fontsize=18)


def viz_binary(
    m, num_dilation_iterations="auto", title=None, ax=None, size_inches=None
):
    """Visualizes a DotPlotMatrix object.

    This is intended to be used on binary DotPlotMatrix objects only.

    Parameters
    ----------
    m: wotplot.DotPlotMatrix
        Matrix object to visualize.

    num_dilation_iterations: int or str
        If this is an integer, then we'll perform this many dilations on the
        matrix before visualization; if this is the string "auto", then we'll
        figure out how many dilations seems right and then do that many.

        For details on dilation, see the OpenCV tutorial linked below. Long
        story short, you can probably leave this at "auto" -- we should only
        need to perform dilations when the input strings get long (more than a
        few hundred nucleotides).

    title: str or None
        If this is not None, then it'll be set as the title of the plot.

    ax: matplotlib.axes.Axes or None
        If this is not None, then we'll add the visualization within this
        Axes object.

    size_inches: (float, float) or None
        If this is not None, then we'll set the resulting visualization's
        Figure to be this large. This will only be used if ax is None: if
        we're drawing this visualization within an existing Axes object, then
        we (as in, this function) don't have control over your plot's size.
        You do ;)

    References
    ----------
    https://docs.opencv.org/4.x/d9/d61/tutorial_py_morphological_ops.html
        OpenCV tutorial on dilation and other morphological operations.

    https://stackoverflow.com/q/44618675
        The answers to this question showed me the light (of image processing
        101) (I didn't know what a dilation was)

    https://matplotlib.org/3.3.4/_modules/matplotlib/axes/_axes.html#Axes.spy
        This function replicates part of what ax.spy() does, with some tweaks.
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
            # sloppy stuff, TODO replace (also check stuff in [500, 2500))
            num_dilation_iterations = 4 + ceil((maxlen - 2500) / 2500)

    if num_dilation_iterations > 0:
        # scipy's dilation function doesn't seem to support sparse matrices...
        # TODO TODO
        matrix_to_show = scipy.sparse.coo_matrix(
            binary_dilation(
                m.mat.toarray(), iterations=num_dilation_iterations
            ).astype(m.mat.dtype)
        )
    else:
        matrix_to_show = m.mat

    return_mplobjs = False
    if ax is None:
        return_mplobjs = True
        fig, ax = pyplot.subplots()
        if size_inches is not None:
            fig.set_size_inches(size_inches)

    _draw_sparse_binary_matrix(ax, matrix_to_show)
    _tidy_dotplot_viz_ax(ax, m)

    if title is not None:
        ax.set_title(title, fontsize=18)

    if return_mplobjs:
        return fig, ax
