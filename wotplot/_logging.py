import time


def get_logger(verbose):
    t0 = time.time()

    def _mlog(s):
        if verbose:
            print(f"{time.time() - t0:,.2f}s: {s}", flush=True)

    return _mlog
