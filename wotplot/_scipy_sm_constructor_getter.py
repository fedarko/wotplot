import scipy


def bail(v):
    raise RuntimeError(
        f'Hey, the SciPy version installed is "{v}", and I don\'t know '
        "how to parse that. Please open a GitHub issue -- sorry!"
    )


def get_vnums(v):
    vnum_parts = v.split(".")
    if len(vnum_parts) >= 2:
        try:
            major_num = int(vnum_parts[0])
            minor_num = int(vnum_parts[1])
        except ValueError:
            bail(v)
        return (major_num, minor_num)
    else:
        bail(v)


def get_sm_constructor():
    """Figures out what sparse matrix constructor to use.

    This is done based on the version of SciPy installed. It will probably
    break eventually, but for now it works and allows us to support Python
    3.6 through 3.11.

    Note about testing: as of writing, if you run wotplot's test suite once,
    then one of the branches in this function will be uncovered (depending
    on what version of SciPy you have installed). However, we use GitHub
    Actions to test wotplot on multiple versions of Python (causing different
    SciPy versions to be installed) -- since the coverage reports from these
    test runs are merged together in CodeCov, this ensures that both branches
    of this function are covered after running the full CI.
    """
    major_num, minor_num = get_vnums(scipy.__version__)
    if major_num >= 2 or (major_num == 1 and minor_num >= 8):
        from scipy.sparse import coo_array

        return coo_array
    else:
        from scipy.sparse import coo_matrix

        return coo_matrix
