from ._make import _make


class DotPlotMatrix:
    """A dot plot matrix showing the relationship between two sequences.

    Attributes
    ----------
    mat (scipy.sparse.coo_array or scipy.sparse.coo_matrix)
        Representation of the dot plot matrix. The exact type will depend on
        the version of SciPy installed.

    k (int)
        k-mer size.

    yorder (str)
        Either "BT" or "TB", indicating the direction that s2 is ordered on the
        y-axis. "BT" indicates "bottom-to-top", and "TB" indicates
        "top-to-bottom".

    grid (bool)
        If True, indicates that the matrix was constructed from collections of
        sequences (rather than from just one sequence on the x-axis and one
        on the y-axis).

    xseqct (int)
        If grid is False, this will be 1. If grid is True, this will be the
        number of sequences in s1.

    yseqct (int)
        If grid is False, this will be 1. If grid is True, this will be the
        number of sequences in s2.

    xlen (int)
        The total length of the sequence(s) along the x-axis.

    ylen (int)
        The total length of the sequence(s) along the y-axis.
    """

    def __init__(
        self,
        s1,
        s2,
        k,
        yorder="BT",
        suff_only=False,
        grid=False,
        verbose=False,
    ):
        """Initializes the dot plot matrix.

        Parameters
        ----------
        s1: str or other str-like object; OR (if grid) collection of these

        s2: str or other str-like object; OR (if grid) collection of these
            The sequences from which we'll create a dot plot. s1 will be on the
            horizontal axis (left to right) and s2 will be on the vertical axis
            (either bottom to top or top to bottom; the order is determined by
            the "yorder" parameter).

            s1 and s2 can be non-str objects (e.g. skbio.DNA), but they both
            must: have length > 0, be convertable to str using str(), and -- in
            their converted-to-str forms -- only contain DNA nucleotides (A, C,
            G, T). We'll always store these sequences as string attributes of
            the newly initialized object.

            If grid is set to True, then we expect s1 and s2 to now be
            collections (lists, tuples, ...) of str or str-like objects.

        k: int
            k-mer size to use when creating the dot plot.

        yorder: str
            Should be either "BT" or "TB". "BT" means that s2 will be ordered
            from bottom to top (like how the Bioinformatics Algorithms textbook
            draws dot plots); "TB" means that s2 will be ordered from top to
            bottom (like how Gepard draws dot plots).

        suff_only: bool
            If True, this uses a different method for identifying shared k-mers
            than the default. The default method (suff_only=False, which uses
            pydivsufsort.common_substrings()) is faster but uses more memory;
            the alternative method (suff_only=True, which uses just
            pydivsufsort's implementation of suffix arrays) is slower but --
            from testing -- uses less memory. See wotplot's README for details.

        grid: bool
            If True, then we will assume that s1 and s2 are collections of
            sequences (rather than individual sequences) -- for example,
            paths of contigs.

        verbose: bool
            If True, logs extra information as this computes the matrix. This
            is useful when working with long sequences and/or for performance
            benchmarking.

        References
        ----------
        Initial implementation based on the Shared k-mers Problem in the
        Bioinformatics Algorithms textbook
        (https://www.bioinformaticsalgorithms.org/) by Compeau & Pevzner.
        """
        self.mat, self.xlen, self.ylen = _make(
            s1, s2, k, yorder, suff_only, verbose, grid=grid
        )
        self.k = k
        self.yorder = yorder
        self.grid = grid
        if self.grid:
            self.xseqct = len(s1)
            self.yseqct = len(s2)
        else:
            self.xseqct = self.yseqct = 1

    def __str__(self):
        if self.yorder == "TB":
            ds = ", top \u2192 bottom"
        else:
            ds = ", bottom \u2192 top"
        if self.grid:
            ds += f", {self.xseqct:,} x {self.yseqct:,} seqs"
        return (
            f"DotPlotMatrix(k = {self.k:,}{ds}): "
            f"{self.mat.shape[0]:,} x {self.mat.shape[1]:,}"
        )

    def __repr__(self):
        # This used to satisfy eval(repr(m)) == m (see
        # https://stackoverflow.com/a/2626364), but our use of a SciPy sparse
        # matrix now means that it doesn't. Also, I removed s1 and s2 in their
        # entirety from this, since when these strings are large it overwhelms
        # the output of this.
        g = f", {self.xseqct}x{self.yseqct} seqs" if self.grid else ""
        return (
            f"DotPlotMatrix(mat={repr(self.mat)}, "
            f'k={self.k}, yorder="{self.yorder}"{g})'
        )
