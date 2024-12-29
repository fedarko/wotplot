import numpy as np
from matplotlib import pyplot
from ._make import FWD, REV, BOTH, MATCH
from ._logging import get_logger

# Colormaps based on Figure 6.20 in Chapter 6 of Bioinformatics Algorithms
# (Compeau & Pevzner), ed. 2
NBCMAP_255 = {
    0: [255, 255, 255],
    FWD: [255, 0, 0],
    REV: [0, 0, 255],
    BOTH: [100, 0, 100],
}
NBCMAP_HEX = {
    0: "#ffffff",
    FWD: "#ff0000",
    REV: "#0000ff",
    BOTH: "#640064",
}

DRAW_ORDER = (FWD, REV, BOTH)


def _get_yarr(yorder):
    if yorder == "BT":
        # -->, which gets turned into an up arrow when we rotate the yax label
        yarr = "\u2192"
    elif yorder == "TB":
        # <--, which gets turned into a down arrow when we rotate the yax label
        yarr = "\u2190"
    else:
        raise ValueError(f"Unrecognized yorder?: {yorder}")
    return yarr


def style_viz_ax(ax, m, s1_name, s2_name, title=None):
    """Adjusts the styling of a matplotlib Axes object to make it look nice.

    Parameters
    ----------
    ax: matplotlib.axes.Axes
    m: wotplot.DotPlotMatrix
    s1_name: str
    s2_name: str
    title: str or None

    Description
    -----------
    Makes the following changes:

    - Hides ticks and tick labels.

    - Adds axis labels formatted like "s1 (x nt)" and "s2 (y nt)", with arrows
      indicating sequence directionality relative to the axis. (The names of s1
      and s2 are controlled by the s1_name and s2_name parameters.)

    - If a title is provided, sets it as the title of the object.
    """
    # Hide the ticks / tick labels.
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.set_xlabel(f"{s1_name} ({len(m.s1):,} nt) \u2192", fontsize=18)
    yarr = _get_yarr(m.yorder)
    ax.set_ylabel(f"{s2_name} ({len(m.s2):,} nt) {yarr}", fontsize=18)

    if title is not None:
        ax.set_title(title, fontsize=18)


def _create_fig_and_ax_if_needed(ax=None):
    if ax is None:
        return pyplot.subplots()
    return None, ax


def viz_spy(
    m,
    markersize=0.5,
    binary=False,
    color="black",
    nbcmap=NBCMAP_HEX,
    draw_order=DRAW_ORDER,
    s1_name="$s_1$",
    s2_name="$s_2$",
    title=None,
    ax=None,
    verbose=False,
    **kwargs,
):
    """Visualizes a DotPlotMatrix object using matplotlib's spy().

    This should be much more performant than viz_imshow(). However, the use of
    a fixed markersize can require manual adjustment based on the size of your
    matrix and desired resolution.

    Parameters
    ----------
    m: wotplot.DotPlotMatrix
        Matrix object to visualize.

    markersize: float
        Size of the markers drawn to represent each match cell. You may want to
        adjust this depending on the size of your matrix.

    binary: bool
        You can set this parameter to True to force the visualization of it as
        a binary matrix (i.e. with forward, reverse-complementary, and
        palindromic match cells all being the same color -- see the "color"
        parameter of this function). For large matrices, this can make the
        visualization step slightly faster.

    color: color
        The color to use for each match cell, if the "binary" parameter is
        True. Must be in a format accepted by matplotlib; see
        https://matplotlib.org/stable/gallery/color/color_demo.html for
        details. (Unlike the colors we use in imshow(), RGB triplets where the
        entries range from 0 to 255 are not allowed here.) If the "binary"
        parameter is False, this is unused.

    nbcmap: dict
        Maps 0, 1, -1, and 2 to colors (of the same possible formats as the
        above "color" parameter). Only used if "binary" is False.

    draw_order: iterable
        If "binary" is False, we will draw the matrix's match cells as distinct
        colors by calling spy() multiple times (once per match type). If your
        markersize is large enough that adjacent cells in the matrix can
        overlap, then the order in which we call spy() will impact which colors
        are drawn "on top" of others in the visualization (with later colors
        ending up on top). You can change the drawing order by adjusting this
        parameter (the default draws forward matches, then
        reverse-complementary matches, then palindromic matches). I don't think
        this should make a big difference in most cases.

    s1_name: str
        Will be used as the name of the sequence on the x-axis.

    s2_name: str
        Will be used as the name of the sequence on the y-axis.

    title: str or None
        If this is not None, then it'll be set as the title of the plot.

    ax: matplotlib.axes.Axes or None
        If this is not None, then we'll add the visualization within this
        Axes object.

    verbose: bool
        If True, prints information about time taken. Useful for performance
        benchmarking.

    **kwargs
        Will be passed to spy().

    Returns
    -------
    fig, ax: (matplotlib.figure.Figure, matplotlib.axes.Axes)
        The figure and axes object created by this function. These are only
        returned if an axes object (the ax parameter above) was not provided;
        if an axes object was provided (i.e. "ax is not None"), then we won't
        return anything.
    """
    _mlog = get_logger(verbose)

    fig, ax = _create_fig_and_ax_if_needed(ax)
    if not binary:
        _mlog("binary is not True, so we'll draw matches in different colors.")
        if len(draw_order) != 3 or set(draw_order) != set(DRAW_ORDER):
            raise ValueError(
                f"draw_order must include exactly 3 elements ({FWD}, {REV}, "
                f"and {BOTH} in any order)."
            )
        # With viz_imshow(), we manually add one RGB triplet per cell
        # (for both match and non-match cells), meaning we don't have to do
        # anything special to set the color of zero cells. In this function,
        # however, we need to explicitly say "okay, set the background color to
        # whatever nbcmap[0] is."
        #
        # spy() doesn't, as far as I can tell, have an easy way to set the
        # background color. We can use set_facecolor(), though -- see
        # https://stackoverflow.com/a/23645437.
        #
        # The most common scenario, I think, will be that the user has default
        # matplotlib styles set (in which the default background color is
        # white) and nbcmap[0] is set to its default (also of white). To save
        # time, I check to see that this is the case -- if so, I don't bother
        # calling set_facecolor(). For all other cases, though, I will call
        # set_facecolor(). (We could try to avoid more redundant uses of this
        # function by converting nbcmap[0] to an RGB triplet to simplify
        # comparison with the output of get_facecolor(), but I'm not sure how
        # much time that approach would even save...)
        if nbcmap[0] != NBCMAP_HEX[0] or ax.get_facecolor() != (1, 1, 1, 1):
            _mlog(f"Setting background color to {nbcmap[0]}...")
            ax.set_facecolor(nbcmap[0])
            _mlog("Done setting background color.")

        # PERF: we could speed this up by only calling spy() for match cell
        # types that actually exist in the matrix (e.g. if there aren't any
        # palindromic matches, skip that spy() call). The most efficient way to
        # do this, I think, would be assigning properties to each matrix during
        # construction that describe the match cell types this matrix has --
        # however, calling spy() with an empty set of values is relatively
        # quick, so this isn't super important.
        for val in draw_order:
            _mlog(f'Visualizing "{val}" cells with spy()...')
            # Filter the matrix to just the cells of a certain match type.
            # https://stackoverflow.com/a/22077616
            # This is somewhat inefficient -- ideally we'd do the filtering in
            # one pass, or somehow make use of the information we already have
            # from matrix construction about where these match cells are.
            ax.spy(
                m.mat.multiply(m.mat == val),
                markersize=markersize,
                color=nbcmap[val],
                **kwargs,
            )
            _mlog(f'Done visualizing "{val}" cells.')
    else:
        _mlog("binary is True; visualizing all match cells with spy()...")
        ax.spy(m.mat, markersize=markersize, color=color, **kwargs)
        _mlog("Done visualizing all match cells.")
    _mlog("Slightly restyling the visualization...")
    style_viz_ax(ax, m, s1_name, s2_name, title)
    _mlog("Done.")
    if fig is not None:
        return fig, ax


def _convert_to_colors(dm, nbcmap):
    """Converts a (dense) ndarray representing a dot plot matrix to colors."""
    # Based on https://stackoverflow.com/a/66821752
    cm = np.ndarray(shape=(dm.shape[0], dm.shape[1], 3))
    for i in range(dm.shape[0]):
        for j in range(dm.shape[1]):
            cm[i][j] = nbcmap[dm[i][j]]
    # Without the .astype("uint8") thing here, I get a warning about
    # "Clipping input data" from imshow().
    # This gets rid of the warning: https://stackoverflow.com/a/58646678
    return cm.astype("uint8")


def viz_imshow(
    m,
    binary=False,
    cmap="gray_r",
    nbcmap=NBCMAP_255,
    s1_name="$s_1$",
    s2_name="$s_2$",
    title=None,
    ax=None,
    verbose=False,
    **kwargs,
):
    """Visualizes a DotPlotMatrix object using matplotlib's imshow().

    IMPORTANT NOTE: This will convert the sparse matrix contained in the
    DotPlotMatrix object to a dense format, in order to make it compatible with
    imshow(). This can require a lot of memory if this dot plot matrix
    describes long sequences -- for large matrices (e.g. where the sequences
    have lengths > ~500 nt), I recommend using viz_spy() instead.

    Parameters
    ----------
    m: wotplot.DotPlotMatrix
        Matrix object to visualize.

    binary: bool
        You can set this parameter to True to force the visualization of it as
        a binary matrix (i.e. with forward, reverse-complementary, and
        palindromic match cells all being the same color). If this is set,
        we'll use the "cmap" parameter to assign colors to the matrix.

    cmap: str
        matplotlib colormap; the default of "gray_r" uses white to represent
        0 (empty cells in the dot plot matrix) and black to represent 1
        (match cells in the dot plot matrix). Only used if "binary" is True.

    nbcmap: dict
        Maps 0, 1, -1, and 2 to colors in RGB triplet format (e.g. red is
        [255, 0, 0]). Only used if "binary" is False.

    s1_name: str
        Will be used as the name of the sequence on the x-axis.

    s2_name: str
        Will be used as the name of the sequence on the y-axis.

    title: str or None
        If this is not None, then it'll be set as the title of the plot.

    ax: matplotlib.axes.Axes or None
        If this is not None, then we will add the visualization within this
        Axes object and not bother creating a new figure and axes object.

    verbose: bool
        If True, prints information about time taken. Useful for performance
        benchmarking.

    **kwargs
        Will be passed to imshow().

    Returns
    -------
    fig, ax: (matplotlib.figure.Figure, matplotlib.axes.Axes)
        The figure and axes object created by this function. These are only
        returned if an axes object (the ax parameter above) was not provided;
        if an axes object was provided (i.e. "ax is not None"), then we won't
        return anything.

    References
    ----------
    The default nbcmap is based on Figure 6.20 in Chapter 6 of Bioinformatics
    Algorithms (Compeau & Pevzner), edition 2.
    """
    _mlog = get_logger(verbose)
    fig, ax = _create_fig_and_ax_if_needed(ax)

    _mlog("Converting the matrix to dense format...")
    dense_mat = m.mat.toarray()
    if not binary:
        _mlog("Converting the dense-format matrix from numbers to colors...")
        dense_mat = _convert_to_colors(dense_mat, nbcmap)
        _mlog("Calling imshow()...")
        ax.imshow(dense_mat, **kwargs)
    else:
        _mlog("Binarizing the dense-format matrix...")
        # I also tried out np.nonzero(), but it looks like the method here
        # (using masks) is faster.
        dense_mat[dense_mat != 0] = MATCH
        _mlog("Calling imshow()...")
        # explicitly set vmin / vmax. This accounts for the corner case where
        # the matrix consists entirely of MATCH cells; in this case, the matrix
        # should be drawn correctly as all full, not all empty (see
        # https://github.com/fedarko/wotplot/issues/19)
        ax.imshow(dense_mat, vmin=0, vmax=MATCH, cmap=cmap, **kwargs)
    _mlog("Slightly restyling the visualization...")
    style_viz_ax(ax, m, s1_name, s2_name, title)
    _mlog("Done.")

    # Only return fig and ax if _create_fig_and_ax_if_needed() created them
    # (i.e. the user did not provide their own "ax" object)
    if fig is not None:
        return fig, ax
