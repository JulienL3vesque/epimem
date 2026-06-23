"""Bundled example data, so the quickstart runs without you supplying your own."""
from importlib import resources

import numpy as np


def example_seasons() -> np.ndarray:
    """R `mem`'s `flucyl` demonstration data, as a (weeks, seasons) array.

    33 surveillance weeks by 8 past influenza seasons of weekly ILI rates (cases per 100,000),
    Castilla y Leon. Use it to run the package end to end before plugging in your own data.
    """
    source = resources.files("epimem").joinpath("data/flucyl.csv")
    with resources.as_file(source) as path:
        return np.loadtxt(path, delimiter=",", skiprows=1)
