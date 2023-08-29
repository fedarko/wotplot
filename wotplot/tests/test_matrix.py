# The array and DotPlotMatrix imports are needed for the
# test_eval_of_repr_equals_orig() test -- if these two terms haven't already
# been loaded, eval() will throw an error.
from wotplot import make


def test_matrix_str_BT_binary():
    d = make("ACC", "CCA", 2)
    assert str(d) == "DotPlotMatrix(k = 2, binary, bottom \u2192 top): 2x2"


def test_matrix_str_TB_binary():
    d = make("ACC", "CCA", 2, yorder="TB")
    assert str(d) == "DotPlotMatrix(k = 2, binary, top \u2192 bottom): 2x2"


def test_matrix_str_TB_notbinary():
    d = make("ACC", "CCA", 2, yorder="TB", binary=False)
    assert str(d) == "DotPlotMatrix(k = 2, top \u2192 bottom): 2x2"
