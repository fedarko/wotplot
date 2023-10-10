import pytest
from wotplot import DotPlotMatrix, FWD, REV, BOTH
from wotplot._viz import viz_spy


def test_viz_spy_bad_draw_order():
    # same example as test_make.test_make_palindrome_not_binary()
    dpm = DotPlotMatrix("AATCGATC", "TATCGAT", 6, binary=False)
    exp_err = (
        f"draw_order must include exactly 3 elements ({FWD}, {REV}, and "
        f"{BOTH} in any order)."
    )
    bad_inputs = [
        (BOTH, BOTH, BOTH),
        (FWD,),
        (FWD, REV, BOTH, 3),
        [],
        (REV, FWD),
    ]
    for bi in bad_inputs:
        with pytest.raises(ValueError) as ei:
            viz_spy(dpm, draw_order=bi)
        assert str(ei.value) == exp_err
