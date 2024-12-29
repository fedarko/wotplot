import time
import random
import wotplot as wp
from matplotlib import pyplot
from memory_profiler import memory_usage


def genseq(n):
    # generates a random DNA sequence of length n
    s = ""
    for i in range(n):
        s += random.choice("ACGT")
    return s


def run(s1, s2, k, suff_only, markersize, ax):
    """Creates and visualizes a dot plot matrix.

    Returns the time taken to run this function.
    """
    print(f"suff_only = {suff_only}...", flush=True)
    t0 = time.time()
    m = wp.DotPlotMatrix(s1, s2, k, suff_only=suff_only, verbose=True)
    t1 = time.time()
    print(f"Matrix construction took {t1 - t0:,.2f} sec.", flush=True)
    wp.viz_spy(
        m, markersize=markersize, title=f"$k$ = {k:,}", verbose=True, ax=ax
    )
    t2 = time.time()
    print(f"Visualization took {t2 - t1:,.2f} sec.", flush=True)
    return t2 - t0


def profile_run(s1, s2, k, suff_only, markersize, ax):
    # returns 2-tuple of (max memory usage, time taken)
    # get memory usage -- https://stackoverflow.com/a/75407972
    return memory_usage(
        (run, (s1, s2, k, suff_only, markersize, ax), {}),
        # get the return value from run(), which is the duration
        retval=True,
        # only run once (otherwise the function could be run multiple times if
        # it is fast enough; see
        # https://github.com/pythonprofilers/memory_profiler/issues/279)
        max_iterations=1,
        # For now, only consider the maximum memory use
        max_usage=True,
    )


def print_sep():
    print("-" * 79, flush=True)


def sim(n, k, markersize=0.5, so_vals=[False, True], fig_size_inches=(10, 7)):
    # Generate dot plot visualizations for two randomly-generated
    # sequences of length n, given a k-mer size k. Profile run(s) of
    # wotplot on these sequences (one for each entry in so_vals).
    # methods, configurable via the suff_only parameter of DotPlotMatrix()).

    print("Generating sequences & prepping before benchmarking...", flush=True)
    s1 = genseq(n)
    s2 = genseq(n)

    if len(so_vals) == 1:
        fig, ax = pyplot.subplots()
        # This is a silly hack but I can't think of a better way to write this
        # atm. granted: i am tired
        if so_vals == [False]:
            axD = ax
        elif so_vals == [True]:
            axS = ax
        else:
            raise ValueError("so_vals has 1 value that isn't a bool? huh.")
    elif len(so_vals) == 2:
        assert len(set(so_vals)) == len(so_vals), "why aren't so_vals unique"
        fig, (axD, axS) = pyplot.subplots(1, 2)
    else:
        raise ValueError("no, stop it, so_vals should have 1 or 2 elements")

    tD = tS = 0
    print_sep()

    if False in so_vals:
        memD, tD = profile_run(s1, s2, k, False, markersize, axD)
        axD.set_title(
            f"common_substrings()\n{tD:,.2f} sec; max mem {memD:,.2f} MiB",
            fontsize=15,
        )
        print_sep()

    if True in so_vals:
        memS, tS = profile_run(s1, s2, k, True, markersize, axS)
        axS.set_title(
            f"suff-only\n{tS:,.2f} sec; max mem {memS:,.2f} MiB", fontsize=15
        )
        print_sep()

    print(f"Total time taken: {tD + tS:,.2f} sec.", flush=True)
    if len(so_vals) == 2:
        suptitley = 0.9
    else:
        suptitley = 1.02
    fig.suptitle(f"$n$ = {n:,}; $k$ = {k:,}", fontsize=22, x=0.51, y=suptitley)
    fig.set_size_inches(fig_size_inches)
