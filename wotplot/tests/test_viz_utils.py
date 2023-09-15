import numpy as np
from wotplot import DotPlotMatrix
from wotplot._viz import _convert_to_colors, NBCMAP_255


def test_convert_to_colors_not_binary():
    # same example as test_make.test_make_palindrome_not_binary()
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    cc = _convert_to_colors(dpm.mat.toarray(), NBCMAP_255)
    assert np.array_equal(
        cc,
        np.array(
            [
                [NBCMAP_255[0], NBCMAP_255[2], NBCMAP_255[0]],
                [NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0]],
            ]
        ),
    )


def test_convert_to_colors_binary():
    # same example as test_make.test_make_simple_default()
    # also, i don't think we should ever call _convert_to_colors() on the
    # output of a binary matrix, but it should work anyway
    dpm = DotPlotMatrix("ACGTC", "AAGTC", 2)
    cc = _convert_to_colors(dpm.mat.toarray(), NBCMAP_255)
    assert np.array_equal(
        cc,
        np.array(
            [
                [NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[1]],
                [NBCMAP_255[1], NBCMAP_255[0], NBCMAP_255[1], NBCMAP_255[0]],
                [NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0]],
                [NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0], NBCMAP_255[0]],
            ]
        ),
    )
