import numpy as np

def ride_hann(d):
    f = 0.5 * (1 - np.cos(2 * np.pi * np.linspace(0, d - 1, d) / (d - 1)))
    return f[:, np.newaxis]  # Convert to a column vector

# # Example usage
# d = 10
# hann_window = ride_hann(d)
# print(hann_window)
