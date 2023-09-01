import time
from pydivsufsort import divsufsort
from ._scipy_sm_constructor_getter import get_sm_constructor

NT2COMP = {"A": "T", "C": "G", "T": "A", "G": "C"}
# Appended to the end of a string when we create its suffix array. "$" occurs
# lexicographically before all of the DNA nucleotides, and including this
# character is helpful when creating suffix arrays -- see Chapter 9 of
# "Bioinformatics Algorithms" for details.
ENDCHAR = "$"

MATCH = 1
FWD = 1
REV = -1
BOTH = 2


def rc(seq):
    out = ""
    for i in range(len(seq) - 1, -1, -1):
        out += NT2COMP[seq[i]]
    return out


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


def _get_shared_kmers(s1, s2, k, s1_sa, s2_sa):
    """Finds the start positions of shared k-mers in two strings.

    Does this using suffix arrays, which makes this more memory-efficient than
    a naive approach.

    Parameters
    ----------
    s1: str
    s2: str
        The strings in which we'll search for shared k-mers. We assume that
        both of these strings have lengths >= k. These strings should NOT
        contain the ENDCHAR character.

    k: int
        k-mer size.

    s1_sa: np.ndarray
        Suffix array for (s1 + ENDCHAR).

    s2_sa: np.ndarray
        Suffix array for (s2 + ENDCHAR).

    Returns
    -------
    matches: list of (int, int)
        Each entry in this list is of the format (p1, p2), and corresponds to
        the presence of a matching k-mer at position p1 in s1 and position p2
        in s2. (Both p1 and p2 are zero-indexed.)

    Example
    -------
    >>> s1 = "ACGTC"
    >>> s2 = "AAGTCAC"
    >>> _get_shared_kmers(
    ...     s1, s2, 2, divsufsort(s1 + ENDCHAR), divsufsort(s2 + ENDCHAR)
    ... )
    [(0, 5), (2, 2), (3, 3)]
    """
    # i and j are indices in the two suffix arrays. we'll set them to 1 in
    # order to skip the first entry in the suffix array (which will always
    # correspond to the suffix of just ENDCHAR, i.e. "$").
    i = 1
    j = 1
    last_kmer_index_1 = len(s1) - k
    last_kmer_index_2 = len(s2) - k
    matches = []
    while i < len(s1_sa) and j < len(s2_sa):
        p1 = s1_sa[i]
        p2 = s2_sa[j]
        if p1 > last_kmer_index_1:
            i += 1
            continue
        if p2 > last_kmer_index_2:
            j += 1
            continue
        # if we've made it here, then we know that both i and j correspond to
        # indices of suffixes in the two strings that each contain at least k
        # characters
        k1 = s1[p1 : p1 + k]
        k2 = s2[p2 : p2 + k]
        if k1 == k2:
            matches.append((p1, p2))
            # "Descend" through s1, identifying all matches between these
            # k-mers and k2. this could probably be made more efficient, or the
            # results here could maybe be reused in future steps (TODO).
            next_p1 = p1 + 1
            while s1[next_p1: next_p1 + k] == k2:
                matches.append((next_p1, p2))
                next_p1 += 1
            # Also "descend" through s2
            next_p2 = p2 + 1
            while s2[next_p2: next_p2 + k] == k1:
                matches.append((p1, next_p2))
                next_p2 += 1
            i += (next_p1 - p1)
            j += (next_p2 - p2)
        else:
            # find lexicographically smaller suffix
            # (We can safely make this comparison using k1 and k2 as proxies
            # for their entire suffix, because we will only end up in this
            # branch if k1 != k2)
            if k1 < k2:
                i += 1
            else:
                j += 1
    return matches


def _make(s1, s2, k, yorder="BT", binary=True, verbose=False):
    """Computes a dot plot matrix.

    Parameters
    ----------
    s1: str or other str-like object
    s2: str or other str-like object
    k: int
    yorder: str
    binary: bool
        See DotPlotMatrix.__init__() for details.

    verbose: bool
        If True, prints a lot of information as this computes the matrix.
        Useful for performance benchmarking.

    Returns
    -------
    (mat, ss1, ss2): (sparse matrix, str, str)
        mat is the main result -- it is the sparse representation of the newly
        created dot plot matrix (the exact type varies depending on the version
        of SciPy installed).

        ss1 and ss2 are versions of s1 and s2, respectively, converted to
        strings.
    """
    t0 = time.time()

    def _mlog(s):
        if verbose:
            print(f"{time.time() - t0:,.2f}s: {s}", flush=True)

    # First, verify that the SciPy version installed is good
    smc = get_sm_constructor()

    # Then validate the inputs
    _mlog("validating inputs...")
    _validate_k(k)
    _validate_yorder(yorder)
    ss1 = _validate_and_stringify_seq(s1, k)
    ss2 = _validate_and_stringify_seq(s2, k)

    # Ok, things seem good. Compute suffix arrays, then find shared k-mers.
    _mlog("computing suffix array for s1...")
    ss1_sa = divsufsort(ss1 + ENDCHAR)

    _mlog("computing suffix array for s2...")
    ss2_sa = divsufsort(ss2 + ENDCHAR)

    _mlog("computing ReverseComplement(s2)...")
    rcs2 = rc(ss2)
    _mlog("computing suffix array for ReverseComplement(s2)...")
    rcs2_sa = divsufsort(rcs2 + ENDCHAR)

    # Find k-mers that are shared between both strings (not considering
    # reverse-complementing)
    _mlog("finding shared k-mers between s1 and s2...")
    fwd_matches = _get_shared_kmers(ss1, ss2, k, ss1_sa, ss2_sa)
    _mlog(f'found {len(fwd_matches):,} such "forward" shared k-mers.')

    _mlog("finding shared k-mers between s1 and ReverseComplement(s2)...")
    rev_matches = _get_shared_kmers(ss1, rcs2, k, ss1_sa, rcs2_sa)
    _mlog(f'found {len(rev_matches):,} such "reverse" shared k-mers.')

    # Convert fwd and rev matches to matrix COO format
    cell2val = {}

    # We could remove the "- k + 1" parts here, but then we'd have empty space
    # for all plots where k > 1 (since you can't have e.g. a 2-mer begin in the
    # final row or column). Interestingly, Figure 6.20 in Bioinformatics
    # Algorithms does include this extra empty space, but we'll omit it here
    mat_shape = (len(ss2) - k + 1, len(ss1) - k + 1)

    def get_row(s2p):
        if yorder == "TB":
            return s2p
        elif yorder == "BT":
            return mat_shape[0] - s2p - 1
        else:
            # should never happen
            raise ValueError(f"Unrecognized yorder: {yorder}")

    def cell_already_fwd(ss1p, ss2p):
        coords = (get_row(ss2p), ss1p)
        # abuse boolean short-circuiting to avoid a KeyError.
        #
        # We coooould make cell2val a defaultdict(int) in order to make this
        # check easier, but then every time we'd try to access a cell that
        # doesn't have a nonzero value assigned yet that cell would get
        # assigned a zero value entry in the defaultdict -- which could
        # unnecessarily increase memory in the sparse matrix (I think "explicit
        # zeroes" take up space). SO ANYWAY using a normal dict avoids this
        # problem
        return coords in cell2val and cell2val[coords] == FWD

    def set_nz_val(val, ss1p, ss2p):
        # val should be FWD or REV. we might change it depending on the matrix
        # type (binary or not) and if another value already exists in this cell
        if binary:
            val_to_use = MATCH
        else:
            if val == REV:
                if cell_already_fwd(ss1p, ss2p):
                    val_to_use = BOTH
                else:
                    val_to_use = REV
            else:
                # We assume that this was called on all FWD matches first, then
                # on all REV matches. If you do this out of order then it'll
                # break the palindrome detection, so don't do that >:(
                val_to_use = FWD
        cell2val[(get_row(ss2p), ss1p)] = val_to_use

    _mlog("converting forward match information to COO format...")
    for m in fwd_matches:
        set_nz_val(FWD, *m)

    _mlog("converting reverse match information to COO format...")
    for m in rev_matches:
        ss2p = len(s2) - m[1] - k
        set_nz_val(REV, m[0], ss2p)

    density = 100 * (len(cell2val) / (mat_shape[0] * mat_shape[1]))
    _mlog(f"{len(cell2val):,} match cell(s); {density:.2f}% density.")
    _mlog("converting to COO format inputs...")
    # Match the input data format expected by SciPy of (vals, (rows, cols)):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html
    mat_vals = []
    mat_rows = []
    mat_cols = []
    for (r, c) in cell2val:
        mat_vals.append(cell2val[(r, c)])
        mat_rows.append(r)
        mat_cols.append(c)

    _mlog("creating sparse matrix...")
    mat = smc((mat_vals, (mat_rows, mat_cols)), shape=mat_shape)
    _mlog("done creating the matrix.")
    return mat, ss1, ss2
