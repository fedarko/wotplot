# The array and DotPlotMatrix imports are needed for the
# test_eval_of_repr_equals_orig() test -- if these two terms haven't already
# been loaded, eval() will throw an error.
from wotplot import make


def test_matrix_str():
    d = make("ACC", "CCA", 2)
    assert str(d) == "DotPlotMatrix(k = 2, binary, bottom \u2192 top): 2x2"
