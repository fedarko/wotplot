# The array and DotPlotMatrix imports are needed for the
# test_eval_of_repr_equals_orig() test -- if these two terms haven't already
# been loaded, eval() will throw an error.
import pytest
from wotplot._scipy_sm_constructor_getter import bail


def test_bail():
    with pytest.raises(RuntimeError) as ei:
        bail("sus version")
    assert str(ei.value) == (
        'Hey, the SciPy version installed is "sus version", and I don\'t know '
        "how to parse that. Please open a GitHub issue -- sorry!"
    )
