"""
=========================================
Using RIDE for speech artifact correction
=========================================

blablablubb
"""

# %%  
# # Loading python modules

# %% 
import pandas as pd
from mne import Epochs, events_from_annotations, set_eeg_reference
from mne.io import read_raw
from pipeline.preprocessing import add_heog_veog, correct_ica

from ride import RideCfg, ride_call, correct_trials

# %%
# # Downloading example data

# %% 

#%% 
# # Read RT data

#%%
log_df = pd.read_csv('/Users/kirstenstark/Documents/research/research_projects/lying_lab_EEG/data/raw/pilots/log/2_memory.txt', 
                     sep='\t', encoding='ISO-8859-1')
log_df = log_df[log_df['bedingung_trigger_erp2'] == 221]

# %% 
# # Preprocessing 

# %%
raw = read_raw("/Users/kirstenstark/Documents/research/research_projects/lying_lab_EEG/data/eeg_VP1/pilot/memory/VP0302.vhdr", preload=True)
raw = raw.set_channel_types({'IO1': 'eog', 'A2': 'misc', 'audio': 'misc', 'pulse': 'misc'})
raw.set_montage('easycap-M1')
raw, _ = set_eeg_reference(raw, ref_channels='average')
raw = add_heog_veog(raw)
#raw, _ = correct_ica(raw, n_components=15)
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
# # Inspecting the RIDE results: Plotting the raw and corrected data

# %%
epochs.average().plot(ylim=dict(eeg=[-25, 25]), titles=dict(eeg="Uncorrected data"), spatial_colors=True)
epochs_corr.average().plot(ylim=dict(eeg=[-25, 25]), titles=dict(eeg="RIDE-corrected data"), spatial_colors=True)

# %% 
# # To Do: VERWEIS AUF PIPELINE 
# # To Do: Fix RIDE
# # To Do: uncomment ICA