#!/usr/bin/env python
# Adapted from https://github.com/fedarko/strainFlye/blob/main/setup.py
# (... adapted from https://github.com/biocore/qurro/blob/master/setup.py)

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

description = "Simple package for creating and visualizing dot plot matrices"

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="wotplot",
    # TODO do something fancy so that this can be stored in __version__
    # (see https://stackoverflow.com/q/458550)
    version="0.1.0",
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
    ],
    setup_requires=[
        "cython",
        "numpy",
    ],
    # Based on how Altair splits up its requirements:
    # https://github.com/altair-viz/altair/blob/master/setup.py
    extras_require={
        "dev": [
            "pytest >= 4.2",
            "pytest-cov >= 2.0",
            "flake8",
            # black 22.10 stopped supporting python 3.6; ideally we'd adjust
            # our CI to use later versions of black on the 3.7 build, but as a
            # hacky solution pinning black is fine
            "black < 22.10",
        ],
        "viz": ["matplotlib", "opencv-python"],
    },
    classifiers=classifiers,
    python_requires=">=3.6",
)
