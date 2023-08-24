import numpy as np
import pytest
from wotplot import make


def test_make_simple():
    dpm = make("ACGTC", "AAGTC", 2)
    assert np.array_equal(
        dpm.mat,
        np.array(
            [
                [0, 0, 0, 1],
                [1, 0, 1, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
            ]
        ),
    )
