import numpy as np


def round_like_matlab(x):
    """Round to nearest integer, like MATLAB's round() function."""

    return np.trunc(x + np.copysign(0.5, x)).astype(int)


def mean_nan(data, dim):

    data = np.nan_to_num(data)
    f = np.mean(data, dim)

    return f
