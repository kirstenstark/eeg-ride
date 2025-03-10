def no_peak(data):
    temp = 0
    for j in range(1, len(data) - 1):
        if data[j - 1] < data[j] and data[j] > data[j + 1]:
            temp = 1
    return temp