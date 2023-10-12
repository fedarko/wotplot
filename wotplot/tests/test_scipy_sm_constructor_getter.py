import pytest
from wotplot._scipy_sm_constructor_getter import bail, get_vnums


def test_bail():
    with pytest.raises(RuntimeError) as ei:
        bail("sus version")
    assert str(ei.value) == (
        'Hey, the SciPy version installed is "sus version", and I don\'t know '
        "how to parse that. Please open a GitHub issue -- sorry!"
    )


def test_get_vnums_good():
    assert get_vnums("1.2") == (1, 2)
    assert get_vnums("2.3.1") == (2, 3)
    assert get_vnums("1.50") == (1, 50)
    assert get_vnums("1234.5678") == (1234, 5678)
    assert get_vnums("100.35.291038123") == (100, 35)


def test_get_vnums_bad():
    for bv in ("a.b", "evil", "1 2 3", ""):
        with pytest.raises(RuntimeError) as ei:
            get_vnums(bv)
        assert str(ei.value) == (
            f'Hey, the SciPy version installed is "{bv}", and I don\'t know '
            "how to parse that. Please open a GitHub issue -- sorry!"
        )
