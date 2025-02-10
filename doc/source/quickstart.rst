Quickstart
==========

Here is one quick example for extracting stimulus- and response-related components from the EEG data of a single participant.

To get started, we need to kinds of input data:

* The preprocessed EEG data, typically generated using the `Epochs` object from `MNE-Python <https://mne.tools/stable/index.html>`_.
  This is a data matrix with the shape ``(n_trials, n_channels, n_timepoints)`` that contains the single trial ERP amplitudes in microvolts.
  Note that for RIDE to work correctly, the data should be cleaned, e.g., by removing eye artifacts, filtering, and rejecting high-amplitude artifacts.
  RIDE needs to be run separately for each experimental condition.
  If you have multiple conditions, please subset the data in an additional step. 

* A list or array of reaction times (RTs) for each trial.
  These are typically extracted from a behavioral log file written by the experiment presentation software.
  The length of the RT array has to match the number of epochs in the EEG data.

Here we assume that the epochs have been stored in MNE's ``.fif`` format and the the RTs are stored in a column of the tab-separated log file.

.. code-block:: python

    import pandas as pd
    from mne import read_epochs

    epochs = read_epochs('sub-01_epo.fif')

    log_df = pd.read_csv('sub-01_logfile.txt', sep='\t')
    rt = log_df['RT']

In the next step, we define the parameters for RIDE.

.. code-block:: python

    from ride import RideCfg

    cfg = RideCfg(comp_name=['s', 'r'],
                  comp_twd=[[0, 600], [-300, 300]],
                  comp_latency=[0, rt],
                  sfreq=500)

Here we've specified:

* The names of the RIDE components that should be extracted.
  At the moment, only one stimulus-related (``'s'``) and one response-related (``'r'``) component are supported.

* The time windows (in milliseconds) in which the RIDE are searched for.
  For the ``'s'`` component, this time window is relative to stimulus onset, while for the ``'r'`` component, it is relative to the reaction time of the trial.
  Note that for the ``'r'`` component, the time window should be long enough to cover any response-related artifact (e.g., speech artifact) entirely.

* The single trial latencies of the RIDE components.
  For the ``'s'`` component, this is always zero, as this is when the stimulus happened relative to the epoch.
  For the ``'r'`` component, this is the list or array of single trial reaction times that we have loaded above.

* The sampling frequency of the EEG data.
  If your input EEG data is a NumPy array, you need to know the sampling frequency from the EEG acquisition.
  If your input EEG data is an ``mne.Epochs`` object, you can check the sampling frequency using `epochs.info.sfreq`.

For additional (optional) input arguments, check the function reference for ``ride.RideCfg``.

Finally, we can run the RIDE procedure using our input EEG data and the configuration:

.. code-block:: python

    from ride import ride_call

    results = ride_call(epochs, cfg)

To inspect the result of the RIDE procedure, you can plot them:


.. code-block:: python

    results.plot()

Here you can see the original ERP (all channels) along with the generated RIDE components ``'s'`` (capturing the stimulus-related activity) and ``'r'`` (capturing the response-related activity, e.g., the speech artifact).

We can use the results from RIDE in a few different ways.
The most typical use case is to remove the ``'r'`` component from the single trial data to remove any response-related artifacts (e.g., speech artifacts).
To do this, the package provides the convenience function ``correct_trials``, taking the fitted RIDE results and the original EEG epochs as inputs:

.. code-block:: python

    from ride import correct_trials

    epochs_corr = correct_trials(results, epochs)

The corrected epochs can then be used for further analysis, e.g., ERP analysis, time-frequency analysis, or source localization. 
