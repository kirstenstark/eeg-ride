import numpy as np
from pandas.api.types import is_list_like
from scipy.signal import correlate


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


def xcov(x, y=None, maxlags=None, scaleopt='none'):
    x = np.asarray(x)
    if y is None:
        y = x
    else:
        y = np.asarray(y)
        
    n = x.shape[0]
    
    if maxlags is None:
        maxlags = n - 1
    else:
        maxlags = int(maxlags)
    
    # Subtract mean
    x = x - np.mean(x)
    y = y - np.mean(y)
    
    # Full cross-correlation
    c = correlate(x, y, mode='full')
    
    lags = np.arange(-n + 1, n)
    
    # Select lags
    mid = len(c) // 2
    start = mid - maxlags
    end = mid + maxlags + 1
    c = c[start:end]
    lags = lags[start:end]
    
    # Scale according to scaleopt
    if scaleopt.lower() == 'biased':
        c = c / n
    elif scaleopt.lower() == 'unbiased':
        c = c / (n - np.abs(lags))
    elif scaleopt.lower() == 'coeff':
        norm_factor = np.std(x) * np.std(y) * n
        c = c / norm_factor
    elif scaleopt.lower() == 'none':
        pass
    else:
        raise ValueError(f"Unknown scaleopt: {scaleopt}")
    
    return c

