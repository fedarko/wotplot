# wotplot

<a href="https://github.com/fedarko/wotplot/actions/workflows/main.yml"><img src="https://github.com/fedarko/wotplot/actions/workflows/main.yml/badge.svg" alt="wotplot CI" /></a>
<a href="https://codecov.io/gh/fedarko/wotplot"><img src="https://codecov.io/gh/fedarko/wotplot/branch/main/graph/badge.svg" alt="Code Coverage" /></a>

wotplot is a small Python library for creating and visualizing
[dot plot matrices](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)).

## Quick examples

### Small dataset

This example is adapted from Figure 6.20 (bottom right) in
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org), edition 2.

```python
import wotplot

# Define our dataset
s1 = "AGCAGGAGATAAACCTGT"
s2 = "AGCAGGTTATCTACCTGT"
k = 3

# Create the matrix (the "binary" parameter means we'll distinguish forward,
# reverse-complementary, and palindromic matching k-mers from each other)
m = wotplot.DotPlotMatrix(s1, s2, k, binary=False)

# Convert the matrix to dense format and visualize it using matplotlib's
# imshow() function (for large matrices where dense representations are
impractical, use viz_spy() instead; see below)
wotplot.viz_imshow(m)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/img/small_example_dotplot.png)

Using the default colorscheme,
:large-red-square: red cells :large-red-square: indicate forward matching
_k_-mers; :large-blue-square: blue cells :large-blue-square: indicate
reverse-complementary matching _k_-mers; and
:large-purple-square: purple cells :large-purple-square: indicate palindromic
matching _k_-mers.

### Larger dataset: comparing two _E. coli_ genomes

Using _E. coli_ K-12 ([from this assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/))
and _E. coli_ O157:H7 ([from this assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008865.2/)).
I removed the two plasmid sequences from the O157:H7 assembly.

```python
import wotplot
from matplotlib import pyplot

# (skipping the part where I loaded the genomes into memory as e1s and e2s...)

# Create the matrix
em = wotplot.DotPlotMatrix(e1s, e2s, 20, verbose=True)

# Visualize the matrix using matplotlib's spy() function
fig, ax = pyplot.subplots()
wotplot.viz_spy(
    em, markersize=0.01, title="Comparison of two $E. coli$ genomes ($k$ = 20)", ax=ax
)
ax.set_xlabel(f"$E. coli$ K-12 substr. MG1655 ({len(e1s)/1e6:.2f} Mbp) \u2192")
ax.set_ylabel(f"$E. coli$ O157:H7 str. Sakai ({len(e2s)/1e6:.2f} Mbp) \u2192")
fig.set_size_inches(8, 8)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/img/ecoli_example_dotplot.png)

## More detailed tutorial

Please see [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Tutorial.ipynb).

## Installation

wotplot supports Python ≥ 3.6. You can install it
and its dependencies using [pip](https://pip.pypa.io):

```
pip install git+https://github.com/fedarko/wotplot.git
```

I'll try to put this on PyPI / conda eventually.

## Performance

### Optimizations made so far

I've tried to make this library reasonably performant. The main optimizations
include:

- We use suffix arrays (courtesy of the lovely
  [`pydivsufsort`](https://github.com/louisabraham/pydivsufsort) library) in
  order to reduce the memory footprint of finding shared _k_-mers.

- We store the dot plot matrix in sparse format (courtesy of
  [SciPy](https://docs.scipy.org/doc/scipy/reference/sparse.html)) in order to
  reduce its memory footprint.

- We support visualizing the dot plot matrix's nonzero values using
  matplotlib's [`spy()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.spy.html)
  function, which (at least for large sequences) is faster and more
  memory-efficient than converting the matrix to a
  dense format and visualizing it with something like
  [`imshow()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.imshow.html).

### That being said...

This library could be made a lot more efficient (I've been documenting ideas in
[issue #2](https://github.com/fedarko/wotplot/issues/2)),
but right now it's good enough for my purposes. Feel free to open an issue / make a pull request
if you'd like to speed it up ;)

### Informal benchmarking

See [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Benchmarking.ipynb)
for some very informal benchmarking results performed on a laptop with ~8 GB of RAM.

Even on this system, the library can handle reasonably large sequences: in the biggest example,
the notebook demonstrates computing the dot plot of two random 100 Mbp sequences
(using _k_ = 20) in 54 minutes and 12.45 seconds.
Dot plots of shorter sequences (e.g. 100 kbp or less) usually take only a few seconds to
compute, at least for reasonably large values of _k_.

## Why does this library exist?

1. This library separates the creation and visualization of dot plot matrices. Other tools that I tried produced pretty visualizations, but didn't give me easy access to the underlying matrix.

2. I wanted something that worked well with [matplotlib](https://matplotlib.org), so that I could create and tile lots of dotplots at once in complicated ways.

## Setting up a development environment

First, fork wotplot. Then you can download a copy of the code from your fork and
install wotplot from this code.

The following commands should work on a Unix system; this assumes that you have
Python ≥ 3.6 and pip installed.

```bash
git clone https://github.com/your-github-username-goes-here/wotplot.git
cd wotplot
pip install -e .[dev]
```

## Contact

Feel free to [open an issue](https://github.com/fedarko/wotplot/issues) if you
have questions, suggestions, comments, or anything else.
