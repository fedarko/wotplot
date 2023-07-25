# wotplot

wotplot is a simple Python 3 package for creating and visualizing
[dot plot matrices](https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)).

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
wotplot.viz(mat)
```

(TODO show screenshot of the above)

This example is borrowed from Figure 6.20 (bottom right) in
[_Bioinformatics Algorithms_](https://www.bioinformaticsalgorithms.org), edition 2.

## Why does this package exist?

1. This package separates the creation of a dot plot matrix from the visualization. Other tools that I tried produced pretty visualizations, but didn't give me easy access to the original matrix.

2. I wanted something that worked well with [matplotlib](https://matplotlib.org), so that I could create and tile lots of dotplots at once in complicated ways.
