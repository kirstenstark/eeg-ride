import numpy as np

def nearest_latency(signal, a):
    """
    Finds the nearest peak latency in a signal based on the given index `a`.

    Parameters:
    signal (array-like): The input signal.
    a (int): The reference index.

    Returns:
    int: The nearest peak index.
    """
    signal = np.array(signal)  # Ensure it's a NumPy array
    temp = np.zeros(len(signal), dtype=int)

    # Identify peaks (local maxima)
    for j in range(1, len(signal) - 1):
        if signal[j + 1] < signal[j] and signal[j] > signal[j - 1]:
            temp[j] = 1

    # If no peaks are found, return 'a' as default
    if np.sum(temp) == 0:
        return a

    # Find the nearest peak to 'a'
    temp1 = np.where(temp == 1)[0]
    temp2 = np.abs(temp1 - a)
    f = temp1[np.argmin(temp2)]
    
    return f

# # Example usage:
# signal = [0, 1, 3, 7, 5, 2, 8, 6, 4, 9, 7]
# a = 5
# nearest_peak = nearest_latency(signal, a)
# print(nearest_peak)
