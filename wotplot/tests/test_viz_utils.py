import pytest
import numpy as np
from wotplot import DotPlotMatrix
from wotplot._viz import _convert_to_colors, _get_yarr, NBCMAP_255


def test_convert_to_colors():
    # same example as test_make.test_make_palindrome_not_binary()
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6)
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


def test_get_yarr():
    assert _get_yarr("BT") == "\u2192"
    assert _get_yarr("TB") == "\u2190"
    with pytest.raises(ValueError) as ei:
        _get_yarr("BB")
    assert str(ei.value) == "Unrecognized yorder?: BB"
