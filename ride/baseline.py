import numpy as np


def baseline(x, a=None):

    x[np.isnan(x)] = 0.0

    if a is None:
        temp = np.mean(x, axis=0)
        f = x - temp[np.newaxis, :]
    else:
        temp = np.mean(x[a, :], axis=0)
        f = x - temp[np.newaxis, :]

    return f
