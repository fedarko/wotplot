# wotplot

<a href="https://github.com/fedarko/wotplot/actions/workflows/main.yml"><img src="https://github.com/fedarko/wotplot/actions/workflows/main.yml/badge.svg" alt="wotplot CI" /></a>
<a href="https://codecov.io/gh/fedarko/wotplot"><img src="https://codecov.io/gh/fedarko/wotplot/branch/main/graph/badge.svg" alt="Code Coverage" /></a>

wotplot is a small Python library for creating and visualizing
[dot plot matrices](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)).

## Quick example

```python
import wotplot

# Define our dataset
s1 = "AGCAGGAGATAAACCTGT"
s2 = "AGCAGGTTATCTACCTGT"
k = 3

# Create the dotplot matrix
m = wotplot.DotPlotMatrix(s1, s2, k)

# Visualize the matrix using matplotlib.imshow()
# (For large matrices, I recommend using viz_spy() instead)
wotplot.viz_imshow(m)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/example_dotplot.png)

This example is adapted from Figure 6.20 (bottom right) in
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org), edition 2.

## More detailed tutorial

Please see [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Tutorial.ipynb).

## Installation

Before installation, wotplot requires Python ≥ 3.6, NumPy, and SciPy.
If you have already installed these dependencies, you can install
wotplot using the following command:

```
pip install git+https://github.com/fedarko/wotplot.git
```

I'll try to put this on PyPI / conda eventually.

## Performance

See [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Benchmarking.ipynb)
for some very informal benchmarking results performed on a laptop with ~8 GB of RAM.

Even on this system, the library can handle reasonably large sequences: in the biggest example,
the notebook demonstrates computing the dot plot of two random 100 Mbp sequences
(using _k_ = 20) in 54 minutes and 12.45 seconds.
Dot plots of shorter sequences (e.g. 100 kbp or less) usually take only a few seconds to
compute, at least for reasonably large values of _k_.

This library [could be made a lot more efficient](https://github.com/fedarko/wotplot/issues/2),
but right now it's good enough for my purposes. Feel free to open an issue / make a pull request
if you'd like to speed it up ;)

## Why does this package exist?

1. This package separates the creation and visualization of dot plot matrices. Other tools that I tried produced pretty visualizations, but didn't give me easy access to the underlying matrix.

2. I wanted something that worked well with [matplotlib](https://matplotlib.org), so that I could create and tile lots of dotplots at once in complicated ways.

## Setting up a development environment

First, fork wotplot. Then you can download a copy of the code from your fork and
install wotplot from this code.

The following commands should work on a Unix system; this assumes that you have
conda installed.

```bash
git clone https://github.com/your-github-username-goes-here/wotplot.git
cd wotplot
conda env create -f environment.yml
conda activate wotplot
pip install -e .[dev]
```

## Contact

Feel free to [open an issue](https://github.com/fedarko/wotplot/issues) if you
have questions, suggestions, comments, or anything else.
