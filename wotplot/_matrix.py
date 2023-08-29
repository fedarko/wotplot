class DotPlotMatrix:
    """Stores a constructed dot plot matrix.

    The main reason this object exists is that I got sick and tired of
    constantly keeping track of the original sequences (and k values,
    etc.) from which I constructed the dot plot matrix. This object doesn't do
    much but store this information in one convenient place.

    Attributes
    ==========
    mat (scipy.sparse.coo_array or scipy.sparse.coo_matrix)
        Representation of the dot plot matrix. The exact type will depend on
        the version of SciPy installed.

    s1 (str)
        Sequence shown on the horizontal axis of the matrix.

    s2 (str)
        Sequence shown on the vertical axis of the matrix.

    k (int)
        k-mer length.

    yorder (str):
        Either "BT" or "TB", indicating the direction that s2 is ordered on the
        y-axis. "BT" indicates "bottom-to-top", and "TB" indicates
        "top-to-bottom".

    binary (bool):
        First, some definitions: a cell in the c-th column and r-th row of mat
        describes the relationship between the k-mer starting at position c in
        s1 and the k-mer starting at position r in s2 (although if yorder ==
        "BT" then the positions in s2 will be flipped). Let's refer to these
        k-mers as k1 and k2, respectively.

        If binary is False, then each cell (c, r) in mat can have one of four
        possible values:

        -  2: k1 == k2, and ReverseComplement(k1) == k2
        -  1: k1 == k2, and ReverseComplement(k1) != k2
        - -1: k1 != k2, and ReverseComplement(k1) == k2
        -  0: k1 != k2, and ReverseComplement(k1) != k2

        If binary is True, then (c, r) can only have two possible values (2, 1,
        and -1 are all all labelled as 1):

        -  1: k1 == k2, and/or ReverseComplement(k1) == k2
        -  0: k1 != k2, and ReverseComplement(k1) != k2
    """

    def __init__(self, mat, s1, s2, k, yorder, binary):
        """Initializes the DotPlotMatrix."""
        self.mat = mat
        self.s1 = s1
        self.s2 = s2
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
