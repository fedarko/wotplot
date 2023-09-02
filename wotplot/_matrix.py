from ._make import _make


class DotPlotMatrix:
    """A dot plot matrix showing the relationship between two sequences.

    Attributes
    ----------
    mat (scipy.sparse.coo_array or scipy.sparse.coo_matrix)
        Representation of the dot plot matrix. The exact type will depend on
        the version of SciPy installed.

    s1 (str)
        Sequence shown on the horizontal axis of the matrix.

    s2 (str)
        Sequence shown on the vertical axis of the matrix.

    k (int)
        k-mer size.

    yorder (str):
        Either "BT" or "TB", indicating the direction that s2 is ordered on the
        y-axis. "BT" indicates "bottom-to-top", and "TB" indicates
        "top-to-bottom".

    binary (bool):
        Whether or not mat is binary. If mat is not binary, then it encodes
        forward, reverse-complementary, and palindromic matches separately; if
        mat is binary, then all matches are represented identically.
    """

    def __init__(self, s1, s2, k, yorder="BT", binary=True, verbose=False):
        """Initializes the dot plot matrix.

        Parameters
        ----------
        s1: str or other str-like object

        s2: str or other str-like object
            The sequences from which we'll create a dot plot. s1 will be on the
            horizontal axis (left to right) and s2 will be on the vertical axis
            (either bottom to top or top to bottom; the order is determined by
            the "yorder" parameter).

            s1 and s2 can be non-str objects (e.g. skbio.DNA), but they both
            must: have length > 0, be convertable to str using str(), and -- in
            their converted-to-str forms -- only contain DNA nucleotides (A, C,
            G, T). We'll always store these sequences as string attributes of
            the newly initialized object.

        k: int
            k-mer size to use when creating the dot plot.

        yorder: str
            Should be either "BT" or "TB". "BT" means that s2 will be ordered
            from bottom to top (like how the Bioinformatics Algorithms textbook
            draws dot plots); "TB" means that s2 will be ordered from top to
            bottom (like how Gepard draws dot plots).

        binary: bool
            If True, then the matrix won't distinguish between forward and
            reverse-complementary matches; if False, it will.

            Some definitions: a cell in the c-th column and r-th row of the
            matrix describes the relationship between the k-mer starting at
            position c in s1 and the k-mer starting at position r in s2
            (although if yorder == "BT" then the positions in s2 will be
            flipped). Let's refer to these k-mers as k1 and k2, respectively.

            If binary is False, then each cell (c, r) in the matrix can have
            one of four possible values:

            -  2: k1 == k2, and ReverseComplement(k1) == k2
            -  1: k1 == k2, and ReverseComplement(k1) != k2
            - -1: k1 != k2, and ReverseComplement(k1) == k2
            -  0: k1 != k2, and ReverseComplement(k1) != k2

            If binary is True, then (c, r) can only have two possible values
            (2, 1, and -1 are all all labelled as 1):

            -  1: k1 == k2, and/or ReverseComplement(k1) == k2
            -  0: k1 != k2, and ReverseComplement(k1) != k2

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
        self.mat, self.s1, self.s2 = _make(s1, s2, k, yorder, binary, verbose)
        self.k = k
        self.yorder = yorder
        self.binary = binary

    def __str__(self):
        bs = ""
        if self.binary:
            bs = ", binary"
        if self.yorder == "TB":
            bs += ", top \u2192 bottom"
        else:
            bs += ", bottom \u2192 top"
        return (
            f"DotPlotMatrix(k = {self.k:,}{bs}): "
            f"{self.mat.shape[0]:,} x {self.mat.shape[1]:,}"
        )

    def __repr__(self):
        # This used to satisfy eval(repr(m)) == m (see
        # https://stackoverflow.com/a/2626364), but our use of a SciPy sparse
        # matrix now means that it doesn't. Also, I removed s1 and s2 in their
        # entirety from this, since when these strings are large it overwhelms
        # the output of this.
        return (
            f"DotPlotMatrix(mat={repr(self.mat)}, "
            f'k={self.k}, yorder="{self.yorder}", binary={self.binary})'
        )
