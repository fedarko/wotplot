import scipy


def bail(sv):
    raise RuntimeError(
        f'Hey, the SciPy version installed is "{sv}", and I don\'t know '
        "how to parse that. Please open a GitHub issue -- sorry!"
    )


def get_sm_constructor():
    """Figures out what sparse matrix constructor to use.

    This is done based on the version of SciPy installed. It will probably
    break eventually, but for now it works and allows us to support Python
    3.6 through 3.11.
    """
    sv = scipy.__version__
    vnum_parts = sv.split(".")

    if len(vnum_parts) >= 2:
        try:
            major_num = int(vnum_parts[0])
            minor_num = int(vnum_parts[1])
        except ValueError:
            bail(sv)
        if major_num >= 2 or (major_num == 1 and minor_num >= 8):
            from scipy.sparse import coo_array

            return coo_array
        else:
            from scipy.sparse import coo_matrix

            return coo_matrix
    else:
        bail(sv)
