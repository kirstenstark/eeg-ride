import numpy as np
from mne import Epochs, EpochsArray

from .helpers import round_like_matlab


def correct_trials(results, data, rt=None):
    """Subtracts the 'r' component from the single trial data, shifted by the RTs."""

    if isinstance(data, (Epochs, EpochsArray)):
        is_epochs = True
        epochs = data.copy()
        ch_types = np.array(epochs.get_channel_types())
        eeg_ixs = np.where(ch_types == 'eeg')[0]
        data = epochs.get_data(picks='eeg').copy()
        data = np.swapaxes(data, 0, 2)

    n_trials = data.shape[2]

    if rt is None:
        rt = results.latencies['r'].copy()
        n_rt = len(rt)
        assert n_trials == n_rt, \
            'When passing different single trial data than the ones used ' + \
            'to perform RIDE, you must also pass the corresponding RTs'
    else:
        rt = rt / results.cfg.re_samp
        rt = round_like_matlab(rt - np.median(rt))
        rt = rt * results.cfg.re_samp

    srate = 1000.0 / results.cfg.samp_interval

    data_corr = data - move3(np.repeat(results.comps['r'][:, :, np.newaxis],
                                       n_trials,
                                       axis=2),
                             np.round(rt*srate/1000).astype(int).flatten())

    if not is_epochs:

        return data_corr

    data_corr = np.swapaxes(data_corr, 0, 2)
    epochs._data[:, eeg_ixs, :] = data_corr

    return epochs


def move3(data, latency):
    """Shifts single trial data ('r' component) by the single trial latency."""

    latency = round_like_matlab(latency)
    d1, d2, d3 = data.shape
    temp = np.zeros((d1, d2, d3), dtype=data.dtype)

    left = latency + 1
    left[left <= 0] = 1
    right = d1 + latency
    right[right > d1] = d1

    left1 = -latency + 1
    left1[left1 <= 0] = 1
    right1 = d1 - latency
    right1[right1 > d1] = d1

    for j in np.where((latency > -d1) & (latency <= d1))[0]:
        temp[left[j] - 1:right[j], :, j] = data[left1[j] - 1:right1[j], :, j]

    return temp
