import numpy as np

from .helpers import round_like_matlab


def median_2d(x):

    # !! Median is probably biased if nan trials are present b/c they are set to zero
    # TODO: simulate this and see if it's a problem
    x = x.copy()
    x[np.isnan(x)] = 0.0  # zeros padding
    s = x.shape
    x = np.sort(x, axis=0)
    y = x.flatten(order='F')[(round_like_matlab(
        s[0]/2) + np.arange(0, s[1])*s[0])-1]

    return y
