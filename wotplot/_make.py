import numpy as np
from collections import defaultdict
from ._matrix import DotPlotMatrix

NT2COMP = {"A": "T", "C": "G", "T": "A", "G": "C"}
    
MATCH = 1
FWD = 1
REV = -1
BOTH = 2

def rc(seq):
    out = ""
    for i in range(len(seq) - 1, -1, -1):
        out += NT2COMP[seq[i]]
    return out


def get_kmer_dd(s, k):
    """Maps each k-mer in a string to a list of start positions of this k-mer."""
    kmers = defaultdict(list)
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
        kmers[kmer].append(i)
    return kmers


def _validate_seq(seq):
    if len(seq) == 0:
        raise ValueError("Input sequence must have length > 0")

def _validate_k(k):
    if k < 1:
        raise ValueError("k must be >= 1")

def _validate_yorder(yorder):
    if yorder not in ("BT", "TB"):
        raise ValueError("yorder must be 'BT' or 'TB'")

def make(s1, s2, k, yorder="BT", binary=True):
    """Computes a dot plot matrix.
    
    Parameters
    ----------
    s1: str
    s2: str
        The strings from which we'll create a dot plot. s1 will be on the
        horizontal axis (left to right) and s2 will be on the vertical axis
        (either bottom to top or top to bottom; the order is determined by
        the "yorder" parameter).
    k: int
        The value of k to use when creating the dot plot.
    yorder: str
        Should be either "BT" or "TB". "BT" means that s2 will be ordered from bottom to top
        (like how the Bioinformatics Algorithms textbook draws dot plots); "TB" means that
        s2 will be ordered from top to bottom (like how Gepard draws dot plots).
    binary: bool
        If True, then the output matrix won't distinguish between forward and
        reverse-complementary matches; if False, it will. See the "Returns"
        section for more detail.

    Returns
    -------
    np.ndarray
        A matrix containing len(s2) rows and len(s1) columns. For a cell in the c-th column
        and r-th row, there are four possible values:

            -  0: No shared k-mers exist that begin at position c in s1 and at position r in s2
            -  1: There is a shared k-mer that begins at position c in s1 and at position r in s2
            - -1: There is a shared k-mer whose RC begins at position c in s1 and at position r in s2
            -  2: There is a palindromic shared k-mer that begins at position c in s1 and at position r in s2

        If "binary" is True, then the output matrix will represent 1,
        -1, and 2 as just 1; if "binary" is False, then the matrix will include
        1, -1, and 2 as distinct values.
    
    References
    ----------
    Based on the Shared k-mers Problem in the Bioinformatics Algorithms
    textbook (https://www.bioinformaticsalgorithms.org/) by Compeau & Pevzner.
    """
    _validate_seq(s1)
    _validate_seq(s2)
    _validate_k(k)
    _validate_yorder(yorder)

    s1_kmers = get_kmer_dd(s1, k)
    s2_kmers = get_kmer_dd(s2, k)

    # We could remove the "- k + 1" parts here, but then we'd have empty space
    # for all plots where k > 1 (since you can't have e.g. a 2-mer begin in the
    # final row or column). Interestingly, Figure 6.20 in Bioinformatics
    # Algorithms does actually include this extra empty space, but I think here
    # it is okay to omit it.
    mat = np.zeros((len(s2) - k + 1, len(s1) - k + 1))
    
    def get_row(s2p):
        if yorder == "TB":
            return s2p
        elif yorder == "BT":
            return mat.shape[0] - s2p - 1
        else:
            # should never happen
            raise ValueError(f"Unrecognized yorder: {yorder}")
    
    # Find k-mers that are shared between both strings (not considering
    # reverse-complementing)
    s1_set = set(s1_kmers.keys())
    s2_set = set(s2_kmers.keys())
    shared_set = s1_set & s2_set
    for shared_kmer in shared_set:
        for s1p in s1_kmers[shared_kmer]:
            for s2p in s2_kmers[shared_kmer]:
                if binary:
                    mat[get_row(s2p)][s1p] = MATCH
                else:
                    mat[get_row(s2p)][s1p] = FWD
                
    # Find k-mers that are shared between both strings, but
    # reverse-complemented
    for s1k in s1_kmers:
        rc_s1k = rc(s1k)
        if rc_s1k in s2_kmers:
            for s1p in s1_kmers[s1k]:
                for s2p in s2_kmers[rc_s1k]:
                    if binary:
                        mat[get_row(s2p)][s1p] = MATCH
                    else:
                        if mat[get_row(s2p)][s1p] == FWD:
                            # If there's both a FWD and RC match here, give it a
                            # unique value
                            mat[get_row(s2p)][s1p] = BOTH
                        else:
                            mat[get_row(s2p)][s1p] = REV

    return DotPlotMatrix(mat, s1, s2, k, binary)
