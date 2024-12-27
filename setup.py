#!/usr/bin/env python
# Adapted from https://github.com/fedarko/strainFlye/blob/main/setup.py
# (... adapted from https://github.com/biocore/qurro/blob/master/setup.py)

import os
from setuptools import find_packages, setup

classifier_str = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Scientific/Engineering :: Visualization
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classifier_str.split("\n") if s]

description = "Small library for creating and visualizing dot plot matrices"

with open("README.md", "r") as f:
    long_description = f.read()

# Extract the version from a special _version.py file. The reason we don't just
# store it in wotplot/__init__.py is that that will break if dependencies
# aren't installed. This is basically option 3 from
# https://packaging.python.org/en/latest/guides/single-sourcing-package-version
# -- also inspired by https://stackoverflow.com/a/7071358
with open(os.path.join("wotplot", "_version.py"), "r") as f:
    # The output of f.read() looks like '__version__ = "0.1.0"\n'.
    # Splitting on "=" gives us ' "0.1.0"\n'; stripping whitespace gives us
    # '"0.1.0"'; and slicing with [1:-1] gives us '0.1.0'.
    version = f.read().split("=")[1].strip()[1:-1]

setup(
    name="wotplot",
    version=version,
    license="BSD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Marcus Fedarko",
    author_email="mfedarko@ucsd.edu",
    maintainer="Marcus Fedarko",
    maintainer_email="mfedarko@ucsd.edu",
    url="https://github.com/fedarko/wotplot",
    packages=find_packages(),
    install_requires=[
        "numpy",
        # earliest version of SciPy I can find that supports python 3 (although
        # the python >= 3.6 requirement will almost certainly result in a
        # higher SciPy version being selected anyway). for reference, this also
        # has scipy.sparse.coo_matrix.
        "scipy >= 0.9",
        "matplotlib",
        "pydivsufsort",
    ],
    # Based on how Altair splits up its requirements:
    # https://github.com/altair-viz/altair/blob/master/setup.py
    extras_require={
        "dev": [
            "pytest >= 4.2",
            "pytest-cov >= 2.0",
            "pytest-mock",
            "flake8",
            # black 22.10 stopped supporting python 3.6; ideally we'd adjust
            # our CI to use later versions of black on the 3.7 build, but as a
            # hacky solution pinning black is fine
            "black < 22.10",
        ],
    },
    classifiers=classifiers,
    python_requires=">=3.6",
)
