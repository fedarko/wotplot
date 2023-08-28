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
    """Maps each k-mer in a string to a list of start positions."""
    kmers = defaultdict(list)
    for i in range(len(s) - k + 1):
        kmer = s[i : i + k]
        kmers[kmer].append(i)
    return kmers


def _validate_and_stringify_seq(seq, k):
    if len(seq) < k:
        raise ValueError(f"Input sequence must have length \u2265 k = {k:,}")
    strseq = str(seq)
    for c in strseq:
        # We cooould get fancy and allow e.g. Uracil, but then that raises the
        # silly question of what if a string contains both T and U??? It's
        # easiest to just mandate that we only take in A/C/G/T strings -- this
        # forces the user to do the conversion, so they can decide what to do.
        if c not in NT2COMP:
            raise ValueError(
                f"Input sequence contains character {c}; only DNA nucleotides "
                "(A, C, G, T) are currently allowed."
            )
    return str(seq)


def _validate_k(k):
    # NOTE: this should be ok even if k is super huge, since python 3 combined
    # "long" and "int" into just "int"; see https://stackoverflow.com/a/2104947
    # We could use "not isinstance(k, int)" if we wanted to be super
    # accommodating, but I don't wanna risk subtle problems since we do pass k
    # over to DotPlotMatrix later
    if k < 1 or type(k) is not int:
        raise ValueError("k must be an integer \u2265 1")


def _validate_yorder(yorder):
    if yorder not in ("BT", "TB"):
        raise ValueError("yorder must be 'BT' or 'TB'")


def make(s1, s2, k, yorder="BT", binary=True):
    """Computes a dot plot matrix.

    Parameters
    ----------
    s1: str or other str-like object
    s2: str or other str-like object
        The sequences from which we'll create a dot plot. s1 will be on the
        horizontal axis (left to right) and s2 will be on the vertical axis
        (either bottom to top or top to bottom; the order is determined by
        the "yorder" parameter). s1 and s2 can be non-str objects (e.g.
        skbio.DNA), but they both must: have length > 0, be convertable to
        str using str(), and -- in their converted-to-str forms -- only contain
        DNA nucleotides (A, C, G, T).
    k: int
        The value of k to use when creating the dot plot.
    yorder: str
        Should be either "BT" or "TB". "BT" means that s2 will be ordered from
        bottom to top (like how the Bioinformatics Algorithms textbook draws
        dot plots); "TB" means that s2 will be ordered from top to bottom (like
        how Gepard draws dot plots).
    binary: bool
        If True, then the output matrix won't distinguish between forward and
        reverse-complementary matches; if False, it will. See the "Returns"
        section for more detail.

    Returns
    -------
    wotplot.DotPlotMatrix

    References
    ----------
    Based on the Shared k-mers Problem in the Bioinformatics Algorithms
    textbook (https://www.bioinformaticsalgorithms.org/) by Compeau & Pevzner.
    """
    _validate_k(k)
    _validate_yorder(yorder)
    ss1 = _validate_and_stringify_seq(s1, k)
    ss2 = _validate_and_stringify_seq(s2, k)

    ss1_kmers = get_kmer_dd(ss1, k)
    ss2_kmers = get_kmer_dd(ss2, k)

    # We could remove the "- k + 1" parts here, but then we'd have empty space
    # for all plots where k > 1 (since you can't have e.g. a 2-mer begin in the
    # final row or column). Interestingly, Figure 6.20 in Bioinformatics
    # Algorithms does actually include this extra empty space, but I think here
    # it is okay to omit it.
    #
    # NOTE: If ss1 and ss2 are both long, this is going to require a horrendous
    # amount of memory. See https://github.com/fedarko/wotplot/issues/2.
    mat = np.zeros((len(ss2) - k + 1, len(ss1) - k + 1))

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
    ss1_set = set(ss1_kmers.keys())
    ss2_set = set(ss2_kmers.keys())
    shared_set = ss1_set & ss2_set
    for shared_kmer in shared_set:
        for ss1p in ss1_kmers[shared_kmer]:
            for ss2p in ss2_kmers[shared_kmer]:
                if binary:
                    mat[get_row(ss2p)][ss1p] = MATCH
                else:
                    mat[get_row(ss2p)][ss1p] = FWD

    # Find k-mers that are shared between both strings, but
    # reverse-complemented
    for ss1k in ss1_kmers:
        rc_ss1k = rc(ss1k)
        if rc_ss1k in ss2_kmers:
            for ss1p in ss1_kmers[ss1k]:
                for ss2p in ss2_kmers[rc_ss1k]:
                    if binary:
                        mat[get_row(ss2p)][ss1p] = MATCH
                    else:
                        if mat[get_row(ss2p)][ss1p] == FWD:
                            # If there's both a FWD and RC match here, give it
                            # a unique value
                            mat[get_row(ss2p)][ss1p] = BOTH
                        else:
                            mat[get_row(ss2p)][ss1p] = REV

    return DotPlotMatrix(mat, ss1, ss2, k, yorder, binary)
