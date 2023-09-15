import numpy as np


def ride_detrend(data, twd):

    # # For debugging
    # data = temp.copy()
    # twd = np.array([51, 110, 288, 348])

    index = np.isnan(data)
    data[index] = 0.0
    d1, d2 = data.shape
    d3 = 1  # MATLAB uses `d1, d2, d3 = data.shape` even though `data` is 2D
    d = [(twd[3] + twd[2]) / 2 - (twd[1] + twd[0]) / 2]

    temp0 = np.tile(np.arange(d1) + 1, [d2, d3]).T

    a = data[np.arange(twd[0], twd[1] + 1), :].mean(axis=0, keepdims=True)
    b = data[np.arange(twd[2], twd[3] + 1), :].mean(axis=0, keepdims=True)
    temp = (b - a) / d
    data = data - temp0 * temp[np.zeros(d1, dtype=int), :]

    temp = data[np.arange(twd[0], twd[1] + 1), :].mean(axis=0, keepdims=True)
    data = data - temp[np.zeros(d1, dtype=int), :]

    f = data.copy()
    f[index] = np.nan

    return f
