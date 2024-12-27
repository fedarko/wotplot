import pytest
from wotplot._make import (
    rc,
    _validate_and_stringify_seq,
    _validate_k,
    _validate_yorder,
    _get_row,
    _fill_match_cells,
)

# technically these are also defined in _make (as of writing), but they are
# intended to be imported from the top-level namespace; so let's do that here.
from wotplot import FWD, REV


def test_rc_good():
    assert rc("ACGTAC") == "GTACGT"
    assert rc("") == ""
    assert rc("A") == "T"
    assert rc("G") == "C"
    assert rc("AT") == "AT"


def test_rc_badchar():
    # in practice, bad characters should already have been caught when make()
    # validates input sequences, and the new method we use for computing
    # reverse-complements (rc()) doesn't disallow non-DNA characters. Let's
    # verify that these "bad characters" don't cause a crash, but keep in mind
    # that these sorts of scenarios should never happen in practice (as of
    # writing).
    assert rc("ACGTM") == "MACGT"
    assert rc("T G") == "C A"


def test_validate_and_stringify_good():
    _validate_and_stringify_seq("A", 1)
    _validate_and_stringify_seq("CTG", 2)


def test_validate_and_stringify_badlen():
    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("", 1)
    assert str(ei.value) == "Input sequence must have length \u2265 k = 1"

    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("CCTGAC", 5678)
    assert str(ei.value) == "Input sequence must have length \u2265 k = 5,678"


def test_validate_and_stringify_badchar():
    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("AC-T", 3)
    assert str(ei.value) == (
        "Input sequence contains character -; only DNA nucleotides (A, C, G, "
        "T) are currently allowed."
    )


def test_validate_k_good():
    for k in (1, 2, 3, 1000, 5000, 10000000):
        _validate_k(k)


def test_validate_k_bad():
    for k in (-1, 0, 1.5, -124324234):
        with pytest.raises(ValueError) as ei:
            _validate_k(k)
        assert str(ei.value) == "k must be an integer \u2265 1"


def test_validate_yorder_good():
    _validate_yorder("BT")
    _validate_yorder("TB")


def test_validate_yorder_bad():
    for b in ("BB", "TT", "", "asdf", 1, 0, None, True, False, "mongus"):
        with pytest.raises(ValueError) as ei:
            _validate_yorder(b)
        assert str(ei.value) == "yorder must be 'BT' or 'TB'"


def test_get_row_good():
    for r in range(0, 5):
        assert _get_row(r, 6, "TB") == r
        assert _get_row(r, 6, "BT") == 5 - r


def test_get_row_bad_yorder():
    with pytest.raises(ValueError) as ei:
        _get_row(3, 6, "BB")
    assert str(ei.value) == "yorder must be 'BT' or 'TB'"


def test_get_row_position_geq_num_rows():
    with pytest.raises(ValueError) as ei:
        _get_row(6, 6, "BT")
    assert str(ei.value) == "s2 pos (6) >= # rows (6)?"

    with pytest.raises(ValueError) as ei:
        _get_row(7, 6, "BT")
    assert str(ei.value) == "s2 pos (7) >= # rows (6)?"

    with pytest.raises(ValueError) as ei:
        _get_row(12345, 6, "BT")
    assert str(ei.value) == "s2 pos (12,345) >= # rows (6)?"


def test_fill_match_cells():
    s1 = "ACGTC"
    s2 = "AAGTCAC"
    md = {}

    _fill_match_cells(s1, s2, 2, md, yorder="TB", binary=False)
    # Note that the coordinates' types will be set based on the types present
    # within the suffix array -- so repr(md) will, for numpy 2, be equal to
    # {(np.int32(5), np.int32(0)): FWD, ...}, at least as of writing.
    # This is not a problem, since np.int32(5) == 5, so the rest of the code
    # should still work.
    assert md == {(5, 0): FWD, (2, 2): FWD, (3, 3): FWD}

    # "Extend" md with reverse-complementary matches
    s2r = rc(s2)
    _fill_match_cells(s1, s2r, 2, md, yorder="TB", binary=False, s2isrc=True)
    assert md == {
        (5, 0): FWD,
        (2, 2): FWD,
        (3, 3): FWD,
        (2, 0): REV,
        (5, 2): REV,
    }


def test_fill_match_cells_redundant_common_substring_rc(mocker):
    gcs = mocker.patch("wotplot._make._get_common_substrings")
    gcs.return_value = [
        (0, 0, 3),
        # Add on an extra unneeded diagonal (the (0, 0, 3) already accounts
        # for matches at (0, 0), (1, 1), and (2, 2), making (1, 1, 2) wholly
        # unnecessary) to verify it doesn't break stuff. (Even though we will
        # process some positions multiple times, we should still keep all of
        # them as REV only -- this is something I was worried about making sure
        # to handle correctly.)
        (1, 1, 2),
    ]

    # I know these two strings have more 1-mer matches than the ones given
    # above (specifically, they have 4 * 3 = 12 distinct 1-mer matches). But,
    # since we've mocked _get_common_substrings() above, there should only be
    # 3 matches now. (Part of the reason I'm adjusting the return value to
    # have fewer matches, besides making the code easier to test, is that this
    # lets us quickly verify that the mocking worked as intended.)
    s1 = "AAAA"
    s2 = "AAA"
    md = {}

    _fill_match_cells(s1, s2, 1, md, yorder="TB", s2isrc=True, binary=False)
    # The positions look a bit different from (0, 0), (1, 1), and (2, 2)
    # because we've set s2isrc=True.
    assert md == {
        (2, 0): REV,
        (1, 1): REV,
        (0, 2): REV,
    }
