import numpy as np
from scipy.fft import fft, ifft


def filtering20(x, a, b):

    f = x
    for j in range(x.shape[1]):
        f[:, j] = filtering10(x[:, j], a, b)

    return f


def filtering10(x, a, b):

    # The MATLAB code does x = x(:) here but this doesn't seem to change anything
    x0 = x.mean()
    x = x - x0
    n = 10 * len(x)
    Y = fft(x, n)

    b = int(b)
    H = np.concatenate([np.zeros((1, a)),
                        np.expand_dims(Y[a: b + 1], axis=0),
                        np.zeros((1, n - b * 2 - 1)),
                        np.expand_dims(Y[n - b: n - a + 1], axis=0),
                        np.zeros((1, a - 1))], axis=1)
    H = H.transpose()

    f = ifft(H, axis=0).real
    f = np.squeeze(f)

    return f[:len(x)]
