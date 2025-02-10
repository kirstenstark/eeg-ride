
# %% [markdown]
# # Example
#
# Here is an example how RIDE can be used at the single trial level to correct speech artifacts in the EEG data.
#
# ## Loading Python modules

# %%
import matplotlib.pyplot as plt
import pandas as pd
from mne import Epochs, events_from_annotations, set_eeg_reference
from mne.io import read_raw
from pipeline.preprocessing import add_heog_veog, correct_ica

from ride import RideCfg, correct_trials, ride_call
from ride.datasets import get_stark2025

# %% [markdown]
# ## Downloading example data
#
# The package comes with a function to download an exemplary dataset of one participant
# from the Stark (2025) study.
# In this EEG experiment, participants performed an object naming task in which a colored frame
# around the object indicated whether participants had to name the shown object (condition
# `simple_naming: 221`) or a previously associated one (condition `memory_naming: 222`).
#
# TThe raw data are stored on the [Open Science Framework](https://osf.io/289ej) and more details
# about the study can be found in [???](https://???).


# %%
log_file, vhdr_file = get_stark2025()

# %% [markdown]
# The paths of the downloaded logfile (`.txt`) and the EEG header file (`.vhdr`)
# can now be fed into the example.

# %%  [markdown]
# ## Read RT data
#
# The log file contains all reaction times (RT) from our single particpant.
# Note that RTs are the latency of the voice key being triggered and are stored in teh column `vkRT`.
#
# Because the RIDE correction is performed for each experimental individually, we here subset to the
# single memory condition (`bedingung_trigger_erp2' == 221`). You can also choose to write a loop for
# all conditions.

# %%
log_df = pd.read_csv(log_file, sep='\t', encoding='ISO-8859-1')
log_df = log_df[log_df['bedingung_trigger_erp2'] == 221]

# %% [markdown]
# ## Preprocessing
#
# Because the RIDE correction is performed on cleaned, preprocessed data, we first perform a basic
# preprocessing steps:
# * Loading the EEG data
# * Setting channel type for non-EEG channels and loading the montage setup for the EEG cap used
#   in the experiment. IO1 is the eye electrode set as EOG. A2, audio, and pulse are additionall
#   channels irrelevant here and therefore set as misc[ellaneous].
# * Re-referencing to common average
# * Performing ICA to correct for eye movement artifacts. We first add virtual VEOG and HEOG channels,
#   perform an independent component analysis using the `fastica` algorithm, and automatically remove
#   components that are likely to reflect blinks or other eye movements. For that, we use the functions
#   `add_heog_veog()` and `correct_ica()` from our general preprocessing pipeline ([Enge, 2021](https://github.com/alexenge/hu-neuro-pipeline),
#   see [Documentation](https://hu-neuro-pipeline.readthedocs.io/en/stable/processing_participant.html)
#   for more details).
# * Applying a bandpass filter between 0.1 and 40 Hz
# * Segmentating to epochs around the relevant stimulus trigger. We here use a time window from -0.2 to 1.5 s,
#   which should include the entire spoken response.
# * Dropping bad epochs and practice trails

# %% tags=["scroll-output"]
raw = read_raw(vhdr_file, preload=True)
raw = raw.set_channel_types(
    {'IO1': 'eog', 'A2': 'misc', 'audio': 'misc', 'pulse': 'misc'})
raw.set_montage('easycap-M1')
raw, _ = set_eeg_reference(raw, ref_channels='average')
raw = add_heog_veog(raw)
raw, _ = correct_ica(raw, n_components=15)
raw = raw.filter(0.1, 40.0)
events, event_id = events_from_annotations(raw)
event_id = {'simple_naming': 221}
epochs = Epochs(raw, events, event_id, tmin=-0.2,
                tmax=1.5, reject=dict(eeg=200e-6))
epochs = epochs[4:]  # remove first 4 epochs as they are practice trials
epochs.metadata = log_df
epochs.drop_bad()
epochs.load_data()

# %% [markdown]
# Feel free to adapt the preprocessing to your needs. However, please note that the
# epochs must be long enough to contain the entire (spoken) reposne for the RIDE
# correction to work properly.

# %% [markdown]
# ## Configuring RIDE
#
# Next, we extract the reponse times `vkRT` from the epochs metadata. These contain
# only the RTs from the good epochs on which we want to base our RIDE correction.
#
# Then, we define the parameters for RIDE:
# * The names of the RIDE components that should be extracted: The stimulus- (`'s'`) and
#   response-related (`'r'`) components.
# * The time windows (in milliseconds) in which the RIDE are searched for. The time window
#   for the `'s'` component is relative to stimulus onset (`comp_twd`: `[0,600]`, `comp_latency`: `0`),
#   while for the `'r'` component, it is relative to the reaction time of the trial
#   (`comp_twd`: `[-400,400]`, `comp_latency`: `rt`). Note that for the `'r'` component, the time window
#   should be long enough to cover any response-related artifact (the speech artifact) entirely.
# * The sampling frequency of the EEG data, extracted from the epochs object.

# %%
rt = list(epochs.metadata['vkRT'])
cfg = RideCfg(comp_name=['s', 'r'],
              comp_twd=[[0, 600], [-400, 400]],
              comp_latency=[0, rt],
              sfreq=epochs.info['sfreq'])

# %% [markdown]
# For additional (optional) input arguments, check the function reference for `ride.RideCfg`.

# %% [markdown]
# ## Running RIDE
#
# Finally, we can run RIDE on the epoched and cleaned EEG data using our configurations.
# The function `ride_call()` returns the RIDE results.

# %% tags=["scroll-output"]
results = ride_call(epochs, cfg)

# %% [markdown]
# ## Inspecting the RIDE results: Plotting the RIDE components
#
# The output of the ride_call function is a dictionary containing, among others, the RIDE components.
# It also contains a plot method to visualize the components.
# The plot shows the an average ERP for all electrodes and the ERP split into the s-component (stimulus-
# related) and the r-component (response-related), i.e. the speech artifact at the median (?) latency.
# Check the ride_call function documentation for further output details.

#  %%
_ = results.plot()

# %% [markdown]
# ## Correcting the original data
#
# The package comes with a separate function to correct the EEG data based on the RIDE results at the
# single trial level. To "clean" the EEG from the speech artifact, the function `correct_trials()`
# subtracts the 'r' component from the single trial data, shifted by the RTs.

# %%
epochs_corr = correct_trials(results, epochs)

# %% [markdown]
# ## Inspecting the RIDE results
#
# ### Plotting the raw and corrected data for one exemplary epoch
#
# Now we can compare the raw and corrected data for one exemplary epoch.

# %%
data = epochs.pick('Fp1')[0].get_data().squeeze()
data_corr = epochs_corr.pick('Fp1')[0].get_data().squeeze()
times = epochs.times

_ = plt.plot(times, data)
_ = plt.plot(times, data_corr)
_ = plt.legend(['uncorrected', 'corrected'])
_ = plt.xlabel('Time (s)')
_ = plt.ylabel('Fp1 Amplitude (V)')

# %% [markdown]
# ## Use RIDE in our pipeline
#
# Instead of manually coding your preprocessing and RIDE correction, you can also use our [hu-neuro-pipeline](https://github.com/alexenge/hu-neuro-pipeline)
# by Enge (2021), based on [Frömer et al, 2018](https://doi.org/10.3389/fnins.2018.00048).
# If you add the argument `ride_perform=True`, this will perform the preprocessing and
# the RIDE correction in one go. Check the [documentation] (https://hu-neuro-pipeline.readthedocs.io/en/stable/inputs_py.html) for managing the further
# pipeline options
