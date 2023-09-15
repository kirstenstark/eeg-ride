import numpy as np


def ride_tukey(n, r):
    n = int(n)
    t = np.linspace(0, 1, n)
    # Define period of the taper as 1/2 period of a sine wave.
    per = r/2
    tl = int(np.floor(per*(n-1))+1)
    th = n-tl+1
    # Window is defined in three sections: taper, constant, taper
    f = np.concatenate([(1+np.cos(np.pi/per*(t[np.arange(tl)] - per)))/2,
                       np.ones((th-tl-1)),
                       (1+np.cos(np.pi/per*(t[th-1:] - 1 + per)))/2])

    if n == 1:
        f = 1

    return f
