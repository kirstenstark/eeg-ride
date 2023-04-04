import matplotlib.pyplot as plt
import numpy as np


def generate_sine_data(
        s_duration=500.0, s_amplitude=15.0, r_duration=200.0, r_amplitude=5.0,
        r_latencies=sorted(np.random.normal(250.0, 50.0, size=100))):
    """
    Generates sine wave data.

    The S component is one long bump (half sine wave) with a fixed onset
    latency of zero. R is a short bump (half sine wave) riding on top of S
    with a variable latency.

    See Figure 1 in Ouyang et al. (2015), https://doi.org/10.1016/j.jneumeth.2014.10.009

    Parameters
    ----------
    s_duration : float
        The duration (in ms) of the S component. This will be the same as the
        duration (first dimension) of the generated data.
    s_amplitude : float
        The amplitude of the S component in (in microvolts).
    r_duration : float
        The duration (in ms) of the R component. Should be shorter than
        `s_duration`.
    r_amplitude : float
        The amplitude of the S component in (in microvolts).
    r_latencies : iterable
        The single trial (peak) latencies of the R component. The length
        of this object will determine the number of trials (third dimension)
        of the generated data.

    Returns
    -------
    data : ndarray (n_times, n_channels, n_trials)
        The generated data (sum of S and single-trial R components). The
        number of channels (second dimension) will always be one.
    """

    times_s = np.arange(s_duration)
    s = s_amplitude * np.sin(np.pi * times_s / s_duration)

    times_r = np.arange(r_duration)
    r = r_amplitude * np.sin(np.pi * times_r / r_duration)

    n_times = int(s_duration)
    n_channels = 1
    n_trials = len(r_latencies)
    data = np.zeros((n_times, n_channels, n_trials))

    channel_ix = 0
    for trial_ix, r_latency in enumerate(r_latencies):
        data[:, channel_ix, trial_ix] = s
        r_start_ix = int(r_latency - r_duration / 2.0)
        r_end_ix = int(r_latency + r_duration / 2.0)
        data[r_start_ix:r_end_ix, channel_ix, trial_ix] += r

    return data, r_latencies


def plot_input_data(data):
    """
    Plots single trial and averaged input data before RIDE.

    TODO: Create one subplot for each channel instead of just plotting the
    first channel.

    Parameters
    ----------
    data : ndarray (n_times, n_channels, n_trials)
        The single trial input data.
    """

    for trial_ix in range(data.shape[2]):
        plt.plot(data[:, 0, trial_ix], color='black', alpha=0.1)

    erp = data.mean(axis=2)
    plt.plot(erp[:, 0], linewidth=4.0, color='red')

    plt.xlabel('Time (ms)')
    plt.ylabel('Amplitude (µV)')
