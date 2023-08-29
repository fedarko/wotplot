import numpy as np
import pytest
from wotplot import DotPlotMatrix


def test_make_simple_default():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2)
    assert np.array_equal(
        dpm.mat.toarray(),
        # The use of a trailing comma on the last line tells black to keep this
        # as a matrix, rather than compress it into a single line. From
        # https://stackoverflow.com/a/66839521.
        np.array(
            [
                [0, 0, 0, 1],
                [1, 0, 1, 0],
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
                [1, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        ),
    )


def test_make_simple_not_binary():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2, binary=False)
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


def test_make_simple_yorder_TB_and_not_binary():
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2, yorder="TB", binary=False)
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


def test_make_palindrome_not_binary():
    # AATCGATC
    # 01234567
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    assert np.array_equal(
        dpm.mat.toarray(),
        np.array(
            [
                [0, 2, 0],
                [0, 0, 0],
            ]
        ),
    )


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
