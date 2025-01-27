"""
=========================================
Using RIDE for speech artifact correction
=========================================

blablablubb
"""

# %%  
# # Loading python modules

# %%
import matplotlib.pyplot as plt 
import pandas as pd
from mne import Epochs, events_from_annotations, set_eeg_reference
from mne.io import read_raw
from pipeline.preprocessing import add_heog_veog, correct_ica

from ride import RideCfg, correct_trials, ride_call
from ride.datasets import get_stark2025

# %%
# # Downloading example data

# %% 
log_file, vhdr_file = get_stark2025()

#%% 
# # Read RT data

#%%
log_df = pd.read_csv(log_file, sep='\t', encoding='ISO-8859-1')
log_df = log_df[log_df['bedingung_trigger_erp2'] == 221]

# %% 
# # Preprocessing 
#
# blablabla 

# %%
raw = read_raw(vhdr_file, preload=True)
raw = raw.set_channel_types({'IO1': 'eog', 'A2': 'misc', 'audio': 'misc', 'pulse': 'misc'})
raw.set_montage('easycap-M1')
raw, _ = set_eeg_reference(raw, ref_channels='average')
raw = add_heog_veog(raw)
raw, _ = correct_ica(raw, n_components=15)
raw = raw.filter(0.1, 40.0)
events, event_id = events_from_annotations(raw)
event_id = {'simple_naming': 221}
epochs = Epochs(raw, events, event_id, tmin=-0.2, tmax=1.5, reject=dict(eeg=200e-6))
epochs = epochs[4:] # remove first 4 epochs as they are practice trials
epochs.metadata = log_df
epochs.drop_bad()
epochs.load_data()



# %%
# # Configuring RIDE

# %%
rt = list(epochs.metadata['vkRT'])
cfg = RideCfg(comp_name=['s', 'r'],
              comp_twd=[[0, 600], [-400, 400]],
              comp_latency=[0, rt],
              sfreq=epochs.info['sfreq'])


# %%
# # Running RIDE

# %%
results = ride_call(epochs, cfg)

# %%
# # Inspecting the RIDE results: Plotting the RIDE components

#  %%
results.plot()

# # Correcting the original data

# %%
epochs_corr = correct_trials(results, epochs)

# %%
# # Inspecting the RIDE results: Plotting the raw and corrected data for one exemplary epoch

# %%
data = epochs.pick('Fp1')[0].get_data().squeeze()
data_corr = epochs_corr.pick('Fp1')[0].get_data().squeeze()
times = epochs.times

plt.plot(times,data)
plt.plot(times,data_corr)
plt.legend(['uncorrected', 'corrected'])
plt.xlabel('Time (s)')
plt.ylabel('Fp1 Amplitude (V)')

# %% 
# # To Do: VERWEIS AUF PIPELINE 
# # To Do: Fix RIDE
# # To Do: uncomment ICA
# %%
