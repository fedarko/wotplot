# wotplot

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

Before installation, wotplot requires a Python version ≥ 3.6; NumPy; and
Cython. If you have already installed these dependencies, you can install
wotplot using the following command:

```
pip install git+https://github.com/fedarko/wotplot.git#egg=wotplot[viz]
```

I'll try to put this on PyPI / conda eventually.

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

First, fork wotplot. Then download the code from your fork.
Something like the following should work; this assumes that you have conda
installed.

```bash
git clone https://github.com/your-github-username-goes-here/wotplot.git
cd wotplot
conda env create -f environment.yml
conda activate wotplot
pip install -e .[dev,viz]
```

## Contact

Feel free to [open an issue](https://github.com/fedarko/wotplot/issues) if you
have any questions, suggestions, comments, or anything else.
