import numpy as np
from scipy.interpolate import CubicSpline


def interp2d(data, x1, x2):

    temp = np.empty((len(x2), data.shape[1]))

    for j in np.arange(data.shape[1]):
        cs = CubicSpline(x1, data[:, j], axis=0)
        temp[:, j] = cs(x2)

    return temp
