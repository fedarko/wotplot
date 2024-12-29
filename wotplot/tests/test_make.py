import pytest
import numpy as np
from wotplot import DotPlotMatrix, MATCH, FWD, REV, BOTH


def test_make_simple_default():
    # Okay, look, MAYBE we could change these if needed. But that's gonna be
    # such a pain. Why would we even do that
    assert FWD == 1
    assert REV == -1
    assert BOTH == 2
    assert MATCH == 1
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2)
    assert np.array_equal(
        dpm.mat.toarray(),
        # The use of a trailing comma on the last line tells black to keep this
        # as a matrix, rather than compress it into a single line. From
        # https://stackoverflow.com/a/66839521.
        np.array(
            [
                [0, 0, 0, 1],
                [-1, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ]
        ),
    )


def test_make_simple_suff_only():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2, suff_only=True)
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 0, 0, 1],
                [-1, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ]
        ),
    )


def test_make_simple_yorder_TB():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2, yorder="TB")
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [-1, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        ),
    )


def test_make_simple_yorder_TB_suff_only():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2, yorder="TB", suff_only=True)
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [-1, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        ),
    )


def test_make_simple_acct():
    # was initially going to be loaded from a FASTA file:
    # https://github.com/fedarko/wotplot/blob/7397b184e8a6dbbc738e1e99f2e2fefd69c9b94f/wotplot/tests/test_make.py
    dpm = DotPlotMatrix("ACCT", "ACCT", 1, yorder="BT")
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [-1, 0, 0, 1],
                [0, 1, 1, 0],
                [0, 1, 1, 0],
                [1, 0, 0, -1],
            ]
        ),
    )


def test_make_simple_acct_ggggg():
    # was initially going to be loaded from a FASTA file:
    # https://github.com/fedarko/wotplot/blob/7397b184e8a6dbbc738e1e99f2e2fefd69c9b94f/wotplot/tests/test_make.py
    dpm = DotPlotMatrix("ACCT", "GGGGG", 1, yorder="BT")
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
            ]
        ),
    )


def test_make_palindrome():
    # AATCGATC
    # 01234567
    assert BOTH == 2
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 2, 0],
                [0, 0, 0],
            ]
        ),
    )


def test_make_palindrome_suff_only():
    # AATCGATC
    # 01234567
    assert BOTH == 2
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, suff_only=True)
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 2, 0],
                [0, 0, 0],
            ]
        ),
    )


def test_make_fancy():
    big = "AGCAGAAAGAGATAAACCTGT"
    dpm = DotPlotMatrix(big, big, 2)
    # This expected output was produced by writing str(m.mat.toarray()) (from
    # an earlier version of wotplot that didn't use the suffix array method,
    # so I had more confidence in it being correct) to a file, and then
    # manipulating this file slightly (removing [] chars, adding \n chars, etc)
    # I think keeping this as a string is nice because it makes it easy to
    # modify this if we need to in the future
    exps = (
        " 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  0  0  1\n"
        " 0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0\n"
        "-1  0  0 -1  0  0  0 -1  0 -1  0  0  0  0  0  0  0  1  0  0\n"
        " 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0\n"
        " 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0 -1\n"
        " 0  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0  0  0\n"
        " 0  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0  0  0\n"
        " 0  0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0\n"
        " 0  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0\n"
        " 0  0  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0\n"
        " 1  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0 -1  0  0\n"
        " 0  0  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0\n"
        " 1  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0 -1  0  0\n"
        " 0  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0  0  0\n"
        " 0  0  0  0  0  1  1  0  0  0  0  0  0  1  1  0  0  0  0  0\n"
        " 0  0  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0\n"
        " 1  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0 -1  0  0\n"
        " 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0\n"
        " 0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
        " 1  0  0  1  0  0  0  1  0  1  0  0  0  0  0  0  0 -1  0  0"
    )
    expsl = exps.split("\n")
    expm = []
    for row in expsl:
        entries = [int(v) for v in row.split()]
        expm.append(entries)
    assert np.array_equal(dpm.mat.toarray(), np.array(expm))


def test_make_bad_chars():
    # bad char in s2
    with pytest.raises(ValueError) as e:
        DotPlotMatrix("ACGTC", "ACDTC", 2)
    assert "Input sequence contains character D" in str(e.value)

    # bad char in s1
    with pytest.raises(ValueError) as e:
        DotPlotMatrix("ACXTC", "ACCTC", 2)
    assert "Input sequence contains character X" in str(e.value)

    # bad char in s1 and s2; s1 should have priority in the error message (but
    # this doesn't really matter that much tbh)
    with pytest.raises(ValueError) as e:
        DotPlotMatrix("ACXTC", "ACDTC", 2)
    assert "Input sequence contains character X" in str(e.value)


def test_make_bad_k():
    for b in (-1, 0, 0.5, 5.4):
        with pytest.raises(ValueError) as e:
            DotPlotMatrix("ACGTC", "ACCTC", b)
        assert str(e.value) == "k must be an integer \u2265 1"


def test_make_suff_only_self_dotplot_reuse_suffixarray(capsys):
    m = DotPlotMatrix(
        "ACCT",
        "ACCT",
        1,
        yorder="BT",
        suff_only=True,
        verbose=True,
    )
    assert np.array_equal(
        m.mat.toarray(),
        np.array(
            [
                [-1, 0, 0, 1],
                [0, 1, 1, 0],
                [0, 1, 1, 0],
                [1, 0, 0, -1],
            ]
        ),
    )
    assert (
        "s1 and s2 are equal, so reusing s1's suffix array for s2..."
        in capsys.readouterr().out
    )


def test_make_suff_only_nonself_dotplot_dont_reuse_suffixarray(capsys):
    m = DotPlotMatrix(
        "ACCT",
        "GGGGG",
        1,
        yorder="BT",
        suff_only=True,
        verbose=True,
    )
    assert np.array_equal(
        m.mat.toarray(),
        np.array(
            [
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
                [0, -1, -1, 0],
            ]
        ),
    )
    assert (
        "s1 and s2 are equal, so reusing s1's suffix array for s2..."
        not in capsys.readouterr().out
    )
