import pydivsufsort
from pydivsufsort import divsufsort
from ._scipy_sm_constructor_getter import get_sm_constructor
from ._logging import get_logger

DNA = "ACGT"
RCDNA = "TGCA"
NT2COMP = str.maketrans(DNA, RCDNA)
# Appended to the end of a string when we create its suffix array. "$" occurs
# lexicographically before all of the DNA nucleotides, and including this
# character is helpful when creating suffix arrays -- see Chapter 9 of
# "Bioinformatics Algorithms" for details.
ENDCHAR = b"$"

MATCH = 1
FWD = 1
REV = -1
BOTH = 2


def rc(seq):
    """Computes the reverse-complement of a DNA string.

    References
    ----------
    Use str.maketrans() to replace nucleotides with their complements; this
    approach (which should be relatively efficient) is based on Devon Ryan's
    suggestion at https://bioinformatics.stackexchange.com/a/3585.
    """
    return seq.translate(NT2COMP)[::-1]


def _validate_and_stringify_seq(seq, k):
    if len(seq) < k:
        raise ValueError(f"Input sequence must have length \u2265 k = {k:,}")
    strseq = str(seq)
    for c in strseq:
        # We cooould get fancy and allow e.g. Uracil, but then that raises the
        # silly question of what if a string contains both T and U??? It's
        # easiest to just mandate that we only take in A/C/G/T strings -- this
        # forces the user to do the conversion, so they can decide what to do.
        if c not in DNA:
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


def _get_suffix_array(seq):
    # We convert the seq to bytes (same as done in pydivsufsort here:
    # https://github.com/louisabraham/pydivsufsort/blob/f4431ee1ea96ee5caf579d9b9e4764636d9cfef1/pydivsufsort/divsufsort.py#L73)
    # in order to prevent a warning about it having to convert the seq to bytes
    # for us. Ideally we wouldn't even work with strings at all (or we'd do
    # this conversion at the start of _make() and use bytes from there on), but
    # I don't really feel like making that change r/n and I don't think it'll
    # make a big difference compared to this tool's other inefficiencies.
    # Although it is still a TODO worth noting (converting both s2 and rc(s2)
    # to bytes separately makes me feel gross).
    return divsufsort(seq.encode("ascii") + ENDCHAR)


def _get_row(position_in_s2, num_rows, yorder):
    """Converts a position in the sequence s2 to a row index in the matrix.

    This is dependent on the yorder of the matrix, and easy to mess up
    accidentally -- which is why I've abstracted this to its own function.

    Parameters
    ----------
    position_in_s2: int
        The position in s2 to convert to a row index. This should be a
        0-indexed position, so it should always be less than num_rows.

    num_rows: int
        The number of rows in the matrix. This should be equal to |s2| - k + 1.

    yorder: str
        Should be "TB" (if s2 is oriented from top -> bottom) or "BT" (if s2 is
        oriented from bottom -> top).

    Returns
    -------
    int
        The row index of position_in_s2 in the matrix.

    Raises
    ------
    ValueError
        - If yorder is not "TB" or "BT".
        - If position_in_s2 >= num_rows.
    """
    _validate_yorder(yorder)
    if position_in_s2 < num_rows:
        if yorder == "TB":
            return position_in_s2
        else:
            return num_rows - position_in_s2 - 1
    else:
        raise ValueError(
            f"s2 pos ({position_in_s2:,}) >= # rows ({num_rows:,})?"
        )


def _fill_match_cells(
    s1, s2, k, cs, md, yorder="BT", binary=True, s2isrc=False
):
    num_rows = len(s2) - k + 1
    for match_run in cs:
        num_cells_matched = match_run[2] - k + 1
        for i in range(num_cells_matched):
            x = match_run[0] + i
            s2p = match_run[1] + i
            if s2isrc:
                s2p = len(s2) - s2p - k
            y = _get_row(s2p, num_rows, yorder)
            pos = (y, x)
            if not binary:
                if s2isrc:
                    if pos in md:
                        md[pos] = BOTH
                    else:
                        md[pos] = REV
                else:
                    md[pos] = FWD
            else:
                md[pos] = MATCH


def _fill_match_cells_old(
    s1, s2, k, s1_sa, s2_sa, md, yorder="BT", binary=True, s2isrc=False
):
    """Finds the start positions of shared k-mers in two strings.

    Does this using suffix arrays, which makes this more memory-efficient than
    a naive approach.

    This used to just return the matching positions, but now its main "output"
    is updating md (a dict that maps matrix position --> match types). This is
    faster than outputting things from here in a nice, easy-to-read format and
    then having to waste time converting that :(

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

    md: dict of (int, int) --> int
        We'll update this dict with keys of the format (p2, p1) (indicating
        that a matching k-mer exists at position p1 in s1 and position p2 in
        s2); these positions correspond to the (row, col) positions of match
        cells in a matrix, taking into account yorder. The values of this dict
        indicate match types (one of {FWD, REV, BOTH, MATCH}).

    yorder: str
        Either "BT" (bottom-to-top) or "TB" (top-to-bottom). See
        DotPlotMatrix.__init__().

    binary: bool
        Either True (only match type used is MATCH) or False (can use FWD, REV,
        or BOTH). See DotPlotMatrix.__init__().

    s2isrc: bool
        Either True (s2 is reverse-complemented) or False (s2 is not
        reverse-complemented).

    Returns
    -------
    None
        (The main "side effect" of this function is updating md; see above.)

    Notes
    -----
    See wotplot/tests/test_make_utils.py for an example of using this function.

    Historical context: I used to have a doctest here, but then numpy 2
    changed how np.int32 types (present in pydivsufsort's output suffix arrays)
    were represented from "7" to "np.int32(7)", which broke the doctest for
    systems running numpy 2. I couldn't find a good general doctest solution
    that worked for both numpy < 2 and numpy 2 (and wasn't ugly).

    Just for reference, this sort of problem w/r/t numpy 2 is documented in
    https://github.com/OSGeo/grass/issues/4100, and the general brittleness
    inherent to doctests is discussed at https://stackoverflow.com/q/13473971.
    """
    # i and j are indices in the two suffix arrays. we'll set them to 1 in
    # order to skip the first entry in the suffix array (which will always
    # correspond to the suffix of just ENDCHAR, i.e. "$").
    i = 1
    j = 1
    last_kmer_index_1 = len(s1) - k
    last_kmer_index_2 = len(s2) - k
    num_rows = len(s2) - k + 1
    while i < len(s1_sa) and j < len(s2_sa):
        p1 = s1_sa[i]
        p2 = s2_sa[j]
        # Since there are only n - k + 1 k-mers in a string of length n, we can
        # ignore the last (k - 1) positions in the string -- there are suffixes
        # that start here, but these suffixes have length < k so we don't care
        # about them.
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
            # "Descend" through s1 and s2, identifying all suffixes where the
            # beginning k-mer matches k1 (and k2, but now we know k1 == k2 so I
            # just use k1 for clarity's sake).
            #
            # NOTE: In the second checks that these while loops make (looking
            # at the beginning k-mer of next_i or next_j), there's the
            # possibility that the positions given by s1_sa[next_i] or
            # s2_sa[next_j] occur after last_kmer_index_1 or last_kmer_index_2,
            # respectively. In these cases, the beginning "k-mer" will be cut
            # off by the end of the string, and will thus by definition be
            # unequal to k1. We could add more explicit checks here, but I'm
            # not sure that the added clarity would be worth the potential
            # slight performance hit. (In big sequences, most positions are ...
            # not near the end of the sequence. damn they should give me a
            # fields medal for that Math Wisdom i just brought into the world.)
            next_i = i + 1
            while (
                next_i < len(s1_sa)
                and s1[s1_sa[next_i] : s1_sa[next_i] + k] == k1
            ):
                next_i += 1

            next_j = j + 1
            while (
                next_j < len(s2_sa)
                and s2[s2_sa[next_j] : s2_sa[next_j] + k] == k1
            ):
                next_j += 1

            # Ok, now we know the "span" of this k-mer in both strings' suffix
            # arrays -- in s1_sa, it's range(i, next_i).
            # (If this k-mer only occurs once in s1, then next_i = i + 1: so
            # list(range(i, next_i)) == [i].)
            #
            # We'll add each match to matches, then jump to just past the ends
            # of these "spans."
            for mi in range(i, next_i):
                x = s1_sa[mi]
                for mj in range(j, next_j):
                    if s2isrc:
                        s2p = len(s2) - s2_sa[mj] - k
                    else:
                        s2p = s2_sa[mj]
                    y = _get_row(s2p, num_rows, yorder)
                    pos = (y, x)
                    if not binary:
                        if s2isrc:
                            if pos in md:
                                md[pos] = BOTH
                            else:
                                md[pos] = REV
                        else:
                            md[pos] = FWD
                    else:
                        md[pos] = MATCH
            i = next_i
            j = next_j
        else:
            # find lexicographically smaller suffix
            # (We can safely make this comparison using k1 and k2 as proxies
            # for their entire suffix, because we will only end up in this
            # branch if k1 != k2)
            if k1 < k2:
                i += 1
            else:
                j += 1


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
    _mlog = get_logger(verbose)

    # First, verify that the SciPy version installed is good
    smc = get_sm_constructor()

    # Then validate the inputs
    _mlog("validating inputs...")
    _validate_k(k)
    _validate_yorder(yorder)
    s1 = _validate_and_stringify_seq(s1, k)
    s2 = _validate_and_stringify_seq(s2, k)

    # Ok, things seem good.

    # We could remove the "- k + 1" parts here, but then we'd have empty space
    # for all plots where k > 1 (since you can't have e.g. a 2-mer begin in the
    # final row or column). Interestingly, Figure 6.20 in Bioinformatics
    # Algorithms does include this extra empty space, but we'll omit it here
    mat_shape = (len(s2) - k + 1, len(s1) - k + 1)

    # TODO: could probably avoid creating "matches" at all, and just populate
    # the matrix directly using the output of common_substrings().
    matches = {}
    _mlog("finding forward matches between s1 and s2...")
    csfwd = pydivsufsort.common_substrings(s1, s2, limit=k)
    _fill_match_cells(s1, s2, k, csfwd, matches, yorder=yorder, binary=binary)
    _mlog(f"found {len(matches):,} forward match cell(s).")
    # I'm not sure if this makes a difference (is Python smart enough to
    # immediately garbage-collect this at this point?), but let's be clear
    del csfwd

    _mlog("computing ReverseComplement(s2)...")
    rcs2 = rc(s2)
    _mlog("finding reverse-complementary matches between s1 and s2...")
    csrev = pydivsufsort.common_substrings(s1, rcs2, limit=k)
    _fill_match_cells(
        s1, rcs2, k, csrev, matches, yorder=yorder, binary=binary, s2isrc=True
    )
    _mlog(f"found {len(matches):,} total match cell(s).")
    del csrev

    density = 100 * (len(matches) / (mat_shape[0] * mat_shape[1]))
    _mlog(f"density = {density:.2f}%.")

    _mlog("converting match information to COO format inputs...")

    # Match the input data format expected by SciPy of (vals, (rows, cols)):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_array.html
    #
    # The reason we don't just output mat_vals, mat_rows, and mat_cols from
    # _find_match_cells() is that we need to be careful about duplicate cells.
    # If binary is False, then we need to do this to identify palindromes; and
    # even if binary is True, we need to do this because including duplicate
    # entries will result in them being summed when creating the matrix
    # (seriously, see the SciPy docs linked above).
    mat_vals = []
    mat_rows = []
    mat_cols = []
    for (r, c) in matches:
        mat_vals.append(matches[(r, c)])
        mat_rows.append(r)
        mat_cols.append(c)

    _mlog("creating sparse matrix from COO format inputs...")
    mat = smc((mat_vals, (mat_rows, mat_cols)), shape=mat_shape)
    _mlog("done creating the matrix.")
    return mat, s1, s2
