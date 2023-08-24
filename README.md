# wotplot

<a href="https://github.com/fedarko/wotplot/actions/workflows/main.yml"><img src="https://github.com/fedarko/wotplot/actions/workflows/main.yml/badge.svg" alt="wotplot CI" /></a>
<a href="https://codecov.io/gh/fedarko/wotplot"><img src="https://codecov.io/gh/fedarko/wotplot/branch/main/graph/badge.svg" alt="Code Coverage" /></a>

wotplot is a simple Python 3 package for creating and visualizing
[dot plot matrices](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)).

It's not very optimized, but it works well enough with relatively small sequences.

## Quick example

```python
import wotplot

# Define our dataset
s1 = "AGCAGGAGATAAACCTGT"
s2 = "AGCAGGTTATCTACCTGT"
k = 3

# Create the dotplot matrix
m = wotplot.make(s1, s2, k)

# Visualize the matrix
wotplot.viz_binary(m)
```

![Output dotplot from the above example](https://github.com/fedarko/wotplot/raw/main/docs/example_dotplot.png)

This example is adapted from Figure 6.20 (bottom right) in
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org), edition 2.

## More detailed tutorial

Please see [this Jupyter Notebook](https://nbviewer.org/github/fedarko/wotplot/blob/main/docs/Tutorial.ipynb).

## Installation

Before installation, wotplot requires a Python version ≥ 3.6 and NumPy.
If you have already installed these dependencies, you can install
wotplot using the following command:

```
pip install git+https://github.com/fedarko/wotplot.git
```

I'll try to put this on PyPI / conda eventually.

### Installation notes

#### Speeding up the installation of `opencv-python`

If you're using Python 3.6, installing `opencv-python` through pip can take a
long time (upwards of 30 minutes).
See [this Stack Overflow thread](https://stackoverflow.com/q/63669752)
for details.

If really you need to use Python 3.6, then I've found that
[specifically running `pip install opencv-python==4.3.0.38`](https://stackoverflow.com/a/63669919)
before you install wotplot fixes this problem. See
[our GitHub Actions workflow](https://github.com/fedarko/wotplot/blob/ce702b63bf790c41d02b0493e3a7eebda6fcec70/.github/workflows/main.yml#L62-L63)
for an example of this.

## Performance

This code is not very optimized, so I don't recommend using it for sequences
longer than a few thousand nucleotides. (On a high-memory system, I've used it
with sequences up to ~50,000 nt; on a laptop, I've used it with sequences up to
~2,000 nt.) See [this issue](https://github.com/fedarko/wotplot/issues/2) for
some discussion of this.

## Why does this package exist?

1. This package separates the creation of a dot plot matrix from the visualization. Other tools that I tried produced pretty visualizations, but didn't give me easy access to the original matrix.

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
