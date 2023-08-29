from wotplot import DotPlotMatrix


def test_matrix_str_BT_binary():
    d = DotPlotMatrix("ACC", "CCA", 2)
    assert str(d) == "DotPlotMatrix(k = 2, binary, bottom \u2192 top): 2 x 2"


def test_matrix_str_TB_binary():
    d = DotPlotMatrix("ACC", "CCA", 2, yorder="TB")
    assert str(d) == "DotPlotMatrix(k = 2, binary, top \u2192 bottom): 2 x 2"


def test_matrix_str_TB_not_binary():
    d = DotPlotMatrix("ACC", "CCA", 2, yorder="TB", binary=False)
    assert str(d) == "DotPlotMatrix(k = 2, top \u2192 bottom): 2 x 2"


def test_matrix_repr():
    d = DotPlotMatrix("ACC", "CCA", 2)
    # I'm not gonna check the entire thing because the repr() of the .mat
    # object will vary depending on its type -- this would make this test super
    # brittle, and although I could account for this I don't think it's worth
    # the trouble
    assert 'k=2, yorder="BT", binary=True)' in repr(d)
