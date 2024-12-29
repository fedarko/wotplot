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


def print_sep():
    print("-" * 79, flush=True)


def sim(n, k, markersize=0.5, fig_size_inches=(10, 7)):
    print("Generating sequences & prepping before benchmarking...")
    # Generate a grid of dot plot visualizations for two randomly-generated
    # sequences of length n, given a k-mer size k.
    s1 = genseq(n)
    s2 = genseq(n)

    fig, (axD, axS) = pyplot.subplots(1, 2)

    print_sep()
    # get memory usage -- https://stackoverflow.com/a/75407972
    memD, tD = memory_usage(
        (run, (s1, s2, k, False, markersize, axD), {}),
        # get the return value from run(), which is the duration
        retval=True,
        # only run once (otherwise the function could be run multiple times if
        # it is fast enough; see
        # https://github.com/pythonprofilers/memory_profiler/issues/279)
        max_iterations=1,
        # For now, only consider the maximum memory use
        max_usage=True,
    )
    print_sep()
    memS, tS = memory_usage(
        (run, (s1, s2, k, True, markersize, axS), {}),
        # get the return value from run(), which is the duration
        retval=True,
        max_iterations=1,
        max_usage=True,
    )

    print_sep()
    axD.set_title(
        f"common_substrings()\n{tD:,.2f} sec; max mem {memD:,.2f} MiB",
        fontsize=18,
    )
    axS.set_title(
        f"suff-only\n{tS:,.2f} sec; max mem {memS:,.2f} MiB", fontsize=18
    )
    print(f"Total time taken: {tD + tS:,.2f} sec.", flush=True)
    fig.suptitle(f"$n$ = {n:,}; $k$ = {k:,}", fontsize=22, x=0.51, y=0.9)
    fig.set_size_inches(fig_size_inches)
