import numpy as np
from matplotlib import pyplot

# Colormap based on Figure 6.20 in Chapter 6 of Bioinformatics Algorithms
# (Compeau & Pevzner), ed. 2
NBCMAP = {
    0: [255, 255, 255],
    1: [255, 0, 0],
    -1: [0, 0, 255],
    2: [100, 0, 100],
}


def style_viz_ax(ax, m, title=None):
    """Adjusts the styling of a matplotlib Axes object to make it look nice.

    Parameters
    ----------
    ax: matplotlib.axes.Axes
    m: wotplot.DotPlotMatrix
    title: str or None

    Description
    -----------
    Makes the following changes:

    - Hides ticks and tick labels.

    - Adds axis labels formatted like "s1 (x nt)" and "s2 (y nt)", with arrows
      indicating sequence directionality relative to the axis.

    - If a title is provided, sets it as the title of the object.
    """
    # Hide the ticks / tick labels.
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

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

    if title is not None:
        ax.set_title(title, fontsize=18)


def _create_fig_and_ax_if_needed(ax=None):
    if ax is None:
        return pyplot.subplots()
    return None, ax


def viz_spy(m, markersize=0.5, color="black", title=None, ax=None, **kwargs):
    """Visualizes a DotPlotMatrix object using matplotlib's spy().

    This should be much more performant than viz_imshow(). However, there are
    some limitations: match cells can only be drawn with a single color, even
    for matrices that are not binary; and you may want to adjust the markersize
    based on the size of your matrix and desired resolution.

    Parameters
    ----------
    m: wotplot.DotPlotMatrix
        Matrix object to visualize.

    markersize: float
        Size of the markers drawn to represent each match cell. You may want to
        adjust this depending on the size of your matrix.

    color: color
        The color to use for each match cell. This should be in a format
        accepted by matplotlib; see
        https://matplotlib.org/stable/gallery/color/color_demo.html for
        details.

    title: str or None
        If this is not None, then it'll be set as the title of the plot.

    ax: matplotlib.axes.Axes or None
        If this is not None, then we'll add the visualization within this
        Axes object.

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

    fig, ax = _create_fig_and_ax_if_needed(ax)
    ax.spy(m.mat, markersize=markersize, color=color, **kwargs)
    style_viz_ax(ax, m, title)
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


def viz_imshow(m, cmap="gray_r", nbcmap=NBCMAP, title=None, ax=None, **kwargs):
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

    cmap: str
        matplotlib colormap; the default of "gray_r" uses white to represent
        0 (empty cells in the dot plot matrix) and black to represent 1
        (match cells in the dot plot matrix). Only used if visualizing a binary
        matrix.

    nbcmap: dict
        Maps 0, 1, -1, and 2 to colors in RGB triplet format (e.g. red is
        [255, 0, 0]). Only used if visualizing a matrix that is not binary.

    title: str or None
        If this is not None, then it'll be set as the title of the plot.

    ax: matplotlib.axes.Axes or None
        If this is not None, then we will add the visualization within this
        Axes object and not bother creating a new figure and axes object.

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
    fig, ax = _create_fig_and_ax_if_needed(ax)

    dense_mat = m.mat.toarray()
    if not m.binary:
        dense_mat = _convert_to_colors(dense_mat, nbcmap)
        ax.imshow(dense_mat, **kwargs)
    else:
        ax.imshow(dense_mat, cmap=cmap, **kwargs)
    style_viz_ax(ax, m, title)

    # Only return fig and ax if _create_fig_and_ax_if_needed() created them
    # (i.e. the user did not provide their own "ax" object)
    if fig is not None:
        return fig, ax
