def find_peak(signal):
    sorted_signal = sorted(enumerate(signal), key=lambda x: x[1])
    return sorted_signal[-1][0]
