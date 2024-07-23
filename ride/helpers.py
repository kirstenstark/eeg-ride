import numpy as np
from pandas.api.types import is_list_like


def round_like_matlab(x):
    """Round to nearest integer, like MATLAB's round() function."""
    if is_list_like(x):
        if sum(np.isnan(x)) > 0:
            raise TypeError('Error: NaN values detected in input array')
    else: 
        if np.isnan(x):
            raise TypeError('Error: NaN values detected in input array')
    return np.trunc(x + np.copysign(0.5, x)).astype(int)


def mean_nan(data, dim):

    data = np.nan_to_num(data)
    f = np.mean(data, dim)

    return f
