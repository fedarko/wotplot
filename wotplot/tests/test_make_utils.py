import pytest
from wotplot._make import (
    rc,
    get_kmer_dd,
    _validate_and_stringify_seq,
    _validate_k,
    _validate_yorder,
)


def test_rc_good():
    assert rc("ACGTAC") == "GTACGT"
    assert rc("") == ""
    assert rc("A") == "T"
    assert rc("G") == "C"
    assert rc("AT") == "AT"


def test_rc_badchar():
    # in practice, bad characters should already have been caught when make()
    # validates input sequences. however, we may as well be careful here...
    with pytest.raises(KeyError):
        rc("ACGTM")


def test_get_kmer_dd_good():
    dd = get_kmer_dd("ACGTTACTC", 2)
    # ACGTTACTC
    # 012345678
    assert dd == {
        "AC": [0, 5],
        "CG": [1],
        "GT": [2],
        "TT": [3],
        "TA": [4],
        "CT": [6],
        "TC": [7],
    }
    assert get_kmer_dd("", 4) == {}
    assert get_kmer_dd("CCT", 3) == {"CCT": [0]}


def test_validate_and_stringify_good():
    _validate_and_stringify_seq("A", 1)
    _validate_and_stringify_seq("CTG", 2)


def test_validate_and_stringify_badlen():
    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("", 1)
    assert str(ei.value) == "Input sequence must have length \u2265 k = 1"

    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("CCTGAC", 5678)
    assert str(ei.value) == "Input sequence must have length \u2265 k = 5,678"


def test_validate_and_stringify_badchar():
    with pytest.raises(ValueError) as ei:
        _validate_and_stringify_seq("AC-T", 3)
    assert str(ei.value) == (
        "Input sequence contains character -; only DNA nucleotides (A, C, G, "
        "T) are currently allowed."
    )


def test_validate_k_good():
    for k in (1, 2, 3, 1000, 5000, 10000000):
        _validate_k(k)


def test_validate_k_bad():
    for k in (-1, 0, 1.5, -124324234):
        with pytest.raises(ValueError) as ei:
            _validate_k(k)
        assert str(ei.value) == "k must be an integer \u2265 1"


def test_validate_yorder_good():
    _validate_yorder("BT")
    _validate_yorder("TB")


def test_validate_yorder_bad():
    for b in ("BB", "TT", "", "asdf", 1, 0, None, True, False, "mongus"):
        with pytest.raises(ValueError) as ei:
            _validate_yorder(b)
        assert str(ei.value) == "yorder must be 'BT' or 'TB'"
