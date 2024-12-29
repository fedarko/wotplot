# wotplot

<a href="https://github.com/fedarko/wotplot/actions/workflows/main.yml"><img src="https://github.com/fedarko/wotplot/actions/workflows/main.yml/badge.svg" alt="wotplot CI" /></a>
<a href="https://codecov.io/gh/fedarko/wotplot"><img src="https://codecov.io/gh/fedarko/wotplot/branch/main/graph/badge.svg" alt="Code Coverage" /></a>
<a href="https://zenodo.org/badge/latestdoi/670484537"><img src="https://zenodo.org/badge/670484537.svg" alt="DOI"></a>
<a href="https://pypi.org/project/wotplot"><img src="https://img.shields.io/pypi/v/wotplot?color=006dad" alt="PyPI" /></a>
<!-- ^ yoinked from pyfastg's github README -->

wotplot is a small Python library for creating and visualizing
[dot plot matrices](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)).

Notably, wotplot creates the _exact_ dot plot matrix, describing
(given some _k_ ≥ 1) every single
[_k_-mer](https://en.wikipedia.org/wiki/K-mer) match between two sequences.
Many tools for visualizing dot plots create only an approximation of this matrix
(containing only the "best" matches) in order to save time; wotplot uses a few
optimizations to make the creation and visualization of the
exact dot plot matrix feasible even for entire prokaryotic genomes. Having this
exact matrix can be useful for a variety of downstream analyses.

## 1. Quick examples

### 1.1. Small dataset

This example is adapted from Figure 6.20 (bottom right) in
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org), edition 2.

```python
import wotplot as wp

# Define our dataset
s1 = "AGCAGGAGATAAACCTGT"
s2 = "AGCAGGTTATCTACCTGT"
k = 3

# Create the matrix
m = wp.DotPlotMatrix(s1, s2, k)

# Convert the matrix to dense format and visualize it using matplotlib's
# imshow() function (for large matrices where dense representations are
# impractical, use viz_spy() instead; see below)
wp.viz_imshow(m)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/img/small_example_dotplot.png)

<!-- Idea of using emojis to represent color c/o https://stackoverflow.com/questions/11509830#comment124410976_41247934 -->
In the default colorscheme
red cells (🟥) indicate forward matches,
blue cells (🟦) indicate reverse-complementary matches,
purple cells (🟪) indicate palindromic matches,
and white cells (⬜) indicate no matches.

### 1.2. Larger dataset: comparing two _E. coli_ genomes

Using _E. coli_ K-12 ([from this assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/))
and _E. coli_ O157:H7 ([from this assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008865.2/)).
I removed the two plasmid sequences from the O157:H7 assembly.

```python
import wotplot as wp

# (skipping the part where I loaded the genomes into memory as e1s and e2s...)

# Create the matrix
# This takes ~35 seconds on a laptop with 8 GB of RAM
em = wp.DotPlotMatrix(e1s, e2s, 20, verbose=True)

# Visualize the matrix using matplotlib's spy() function
# This takes ~1 second on a laptop with 8 GB of RAM
fig, ax = wp.viz_spy(em, markersize=0.01, title="Comparison of two $E. coli$ genomes ($k$ = 20)")
ax.set_xlabel(f"$E. coli$ K-12 substr. MG1655 ({len(e1s)/1e6:.2f} Mbp) \u2192")
ax.set_ylabel(f"$E. coli$ O157:H7 str. Sakai ({len(e2s)/1e6:.2f} Mbp) \u2192")
fig.set_size_inches(8, 8)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/img/ecoli_example_dotplot.png)

## 2. More detailed tutorial

Please see [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Tutorial.ipynb).

## 3. Installation

wotplot supports Python ≥ 3.6. You can install it and its dependencies using
[pip](https://pip.pypa.io):

```bash
pip install wotplot
```

## 4. Performance

### 4.1. Optimizations made so far

I've tried to make this library reasonably performant. The main optimizations
include:

- We use the [`pydivsufsort`](https://github.com/louisabraham/pydivsufsort)
  library to find shared _k_-mers (corresponding to filled-in cells in the
  dot plot matrix). See below for details.

- We store the dot plot matrix in sparse format (courtesy of
  [SciPy](https://docs.scipy.org/doc/scipy/reference/sparse.html)) in order to
  reduce its memory footprint.

- We support visualizing the dot plot matrix's nonzero values using
  matplotlib's [`spy()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.spy.html)
  function, which (at least for large sequences) is faster and more
  memory-efficient than converting the matrix to a
  dense format and visualizing it with something like
  [`imshow()`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.imshow.html).

### 4.2. Two methods for finding shared _k_-mers

In order to create an exact dot plot comparing two sequences, we need to
identify all pairs of shared _k_-mers between these sequences and their
locations within these sequences.
Doing this in a time- and space-efficient way is tricky.

wotplot supports two methods for finding shared _k_-mers, both of which rely on
the excellent [`pydivsufsort`](https://github.com/louisabraham/pydivsufsort) library:

1. **Default**: just uses **[`pydivsufsort.common_substrings()`](https://github.com/louisabraham/pydivsufsort/issues/42)** (faster, but requires more memory)

2. **"suff-only":** uses **`pydivsufsort.divsufsort()`** to compute suffix arrays, then iterates through them (slower, but requires less memory)

If you want to use the suff-only method, you can set `suff_only=True` when creating a
`DotPlotMatrix` object.

#### 4.2.1. Details about these methods

So, I implemented the second method ("suff-only") first, way back in 2023. This method
uses `pydivsufsort.divsufsort()` (which in turn uses
[`libdivsufsort`](https://github.com/y-256/libdivsufsort)) to compute suffix arrays for
each of the input strings, then iterates through these suffix arrays simultaneously to
identify shared _k_-mers. The iterating-through-suffix-arrays part of this method (aka
the stuff that I wrote) is relatively simple, and could probably be made much more efficient.

More recently, I switched to using
[`pydivsufsort.common_substrings()`](https://github.com/louisabraham/pydivsufsort/issues/42)
as the default method for identifying shared _k_-mers. It's much faster than the "suff-only"
method -- for example, computing the dot plot of the _E. coli_ example shown above takes
~3 minutes with the suff-only method, but ~35 seconds with the `common_substrings()` method!

Howeever, I've noticed from benchmarking (see below) that the `common_substrings()` method
has higher memory requirements for long sequences than the suff-only method. For example,
on my laptop with 8 GB of RAM, wotplot using the `common_substrings()` method sometimes crashes when
creating dot plots of two 20 Mbp sequences (although, when it doesn't crash, it can create such
a dot plot in about 62 seconds). The peak memory usage from such a successful run
is ~5,823.75 MiB (aka ~6.11 GB).

wotplot using the suff-only method, for comparison, can create a dot plot of two
150 Mbp (!) sequences on the same laptop -- with peak memory usage of ~2,318.79 MiB
(aka ~2.43 GB). The downside is that it is slow; creating this massive dot plot takes over an hour.

Anyway -- for most use cases, I think `common_substrings()` will be a better choice, so
I'm leaving it as the default. However, if you are working with long sequences on
low-memory systems, you may need to use the suff-only method.

#### 4.2.2. When should I use one method or another?

It depends on how much memory your system has and how long your sequences are. Speaking
very generally, assuming you are on a system with ~8 GB RAM, the default method should be
fine up until your sequences are ~5-10 Mbp each. At that point, if you are okay with taking a
longer time to create your dot plots, you may want to use the suff-only method (or at least
think about it, in case you eventually start running out of memory).

(If you need to create dot plots of (i) very long sequences (ii) on a low-memory system and (iii) you
need to do it as quickly as possible, this library might not be ideal -- since it is creating the
exact dot plot matrix. Using a tool that creates a less granular dot plot might better meet your needs.)

### 4.3. Informal benchmarking

See [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/tree/main/docs/Benchmarking.ipynb).

### 4.4. Some ideas for improving performance

I've documented a few ideas in [issue #2](https://github.com/fedarko/wotplot/issues/2).
Right now this library is fast enough for my purposes, but feel free to open an
issue / make a pull request if you'd like to speed it up :)

## 5. Okay but like why does this library exist

1. This library separates the creation and visualization of dot plot matrices. Other tools that I tried produced pretty visualizations, but didn't give me easy access to the underlying matrix.

2. I wanted something that worked well with [matplotlib](https://matplotlib.org), so that I could create and tile lots of dotplots at once in complicated ways.

## 6. Limitations

- **Performance:** Although I've tried to optimize this tool (see above),
  it definitely isn't the fastest or the most memory-efficient way to
  visualize dot plots. The two obvious reasons for this are that (1) this
  is written in Python, and (2) this is creating the exact dot plot matrix
  rather than a subset of it.

- **Only static visualizations:** The visualization methods included in the
  tool only support the creation of static plots. There are
  [ways to make matplotlib visualizations interactive](https://matplotlib.org/stable/users/explain/interactive.html) (e.g. using
  [`%matplotlib notebook`](https://stackoverflow.com/a/41125787) within a
  Jupyter Notebook, or with the `pyplot.show()` GUI), but (1) I don't
  know enough about these methods to "officially" support them and (2) these
  visualizations will still probably pale in comparison to the outputs of
  dedicated interactive visualization software
  (e.g. [ModDotPlot](https://github.com/marbl/ModDotPlot)).

## 7. Setting up a development environment

First, fork wotplot -- this will make it easy to submit a pull request later.

After you've forked wotplot, you can download a copy of the code from your
fork and install wotplot from this downloaded code. The following commands
should do this; note that these commands assume (1) that you're using a
Unix system and (2) that you have Python ≥ 3.6 and pip installed.

```bash
git clone https://github.com/your-github-username-goes-here/wotplot.git
cd wotplot
pip install -e .[dev]
```

After the above commands, you can check that wotplot was installed successfully
by running its test suite:

```bash
make test
```

## 8. Acknowledgements

The small example given above, and my initial implementation of an algorithm
for computing dot plots, were based on Chapter 6 of
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org)
(Compeau & Pevzner).

The idea of using suffix arrays to speed up dot plot computation is not new; it
is also used in [Gepard](https://cube.univie.ac.at/gepard)
([Krumsiek _et al._ 2007](https://academic.oup.com/bioinformatics/article/23/8/1026/198110)).
(Eventually I moved from directly using suffix arrays to just using the
`pydivsufsort.common_substrings()` algorithm, at least as the default algorithm; although
that is still [using a suffix array under the hood](https://github.com/louisabraham/pydivsufsort/blob/2869020c26022e0f88592e85cdc480856e9856d5/pydivsufsort/wonderstring.py#L128-L157) :)

### 8.1. Dependencies

- [NumPy](https://numpy.org)
- [SciPy](https://scipy.org)
- [`pydivsufsort`](https://github.com/louisabraham/pydivsufsort)
- [matplotlib](https://matplotlib.org)

### 8.3. Testing dependencies

- [pytest](https://docs.pytest.org)
- [pytest-cov](https://github.com/pytest-dev/pytest-cov)
- [pytest-mock](https://github.com/pytest-dev/pytest-mock)
- [flake8](https://flake8.pycqa.org)
- [black](https://github.com/psf/black)

## 9. Contact

Feel free to [open an issue](https://github.com/fedarko/wotplot/issues) if you
have questions, suggestions, comments, etc.
