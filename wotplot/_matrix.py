class DotPlotMatrix:
    """Stores a constructed dot plot matrix.

    The main reason this object exists is that I got sick and tired of
    constantly keeping track of the original sequences (and k values,
    etc.) from which I constructed the dot plot matrix. This object doesn't do
    much but store this information in one convenient place.
    """

    def __init__(self, mat, s1, s2, k):
        self.mat = mat
        self.s1 = s1
        self.s2 = s2
        self.k = k

    def __str__(self):
        return (
            f"DotPlotMatrix(k = {self.k:,}): "
            f"{self.mat.shape[0]:,} x {self.mat.shape[1]:,}"
        )

    def __repr__(self):
        # This satisfies eval(repr(m)) == m; see
        # https://stackoverflow.com/a/2626364.
        return (
            f"DotPlotMatrix(mat={repr(self.mat)}, "
            f's1="{self.s1}", s2="{self.s2}", k={self.k})'
        )
