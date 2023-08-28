import pytest
from collections import defaultdict
from wotplot._make import rc, get_kmer_dd


def test_rc_simple():
    assert rc("ACGTAC") == "GTACGT"
    assert rc("") == ""
    assert rc("A") == "T"
    assert rc("G") == "C"
    assert rc("AT") == "AT"


def test_rc_badchar():
    # in practice, bad characters should already have been caught when make()
    # validates input sequences. however, we may as well be careful here...
    with pytest.raises(KeyError) as ei:
        rc("ACGTM")


def test_get_kmer_dd_simple():
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
