# The array and DotPlotMatrix imports are needed for the
# test_eval_of_repr_equals_orig() test -- if these two terms haven't already
# been loaded, eval() will throw an error.
# We can tell flake8 "it's cool, i really am using these" by using the noqa
# syntax; see https://stackoverflow.com/a/59167435.
from numpy import array  # noqa: F401
from wotplot import DotPlotMatrix, make  # noqa: F401


def test_eval_of_repr_equals_orig():
    dpm = make("ACC", "CCA", 2)
    assert eval(repr(dpm)) == dpm
