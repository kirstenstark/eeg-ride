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





# # Import Python module
# pipeline <- reticulate::import("pipeline")
# 

# # Run the group level pipeline
# res_erp2 <- pipeline$group_pipeline(
#     vhdr_files=EEG_raws,
#     log_files=here::here("data", "raw", "log","with_practice_for_EEG_memory"),
#     output_dir=here::here("data", "eeg_VP1", "export", "memory_ERP2"),
#     # report_dir = here::here("data", "eeg_VP1", "export",
#     #                         "qc_reports_card_ERP1"), # only first time
#     downsample_sfreq = NULL, # default
#    veog_channels = "auto", #c("Fp1", "IO1"), # default
#    heog_channels = "auto", # c("F9", "F10"), # default
#    montage = "easycap-M1", # default
#    bad_channels = list(
#      "VP0301" = c("T7"),
#      # "VP0312" = c("TP9"), # nur wenn zu wenige Segmente
#      # "VP0315" = c("TP9"), # nur wenn zu wenige Segmente
#      # "VP0327" = c("TP9"), # nur wenn zu wenige Segmente
#     "VP0330" = c("TP10")
#      # "VP0332" = c("FC6"), # nur wenn zu wenige Segmente
#      ),
#    besa_files = NULL,
#    ica_method =  "fastica",
#    ica_n_components = 0.99, # prop. of var. explained by ICA comp.
#    highpass_freq =  0.1,
#    lowpass_freq = 40.0,
#     triggers = c(221, 222),
#     triggers_column = "bedingung_trigger_erp2",
#     epochs_tmin = -0.2,
#     epochs_tmax = 2.0,
#     baseline = c(-0.2, 0.0),
#    reject_peak_to_peak = 200.0,
#     components = list(
#         "name" = list("P1", "N1", "N2", "P300", "posP2"),
#         "tmin" = list(0.080, 0.130, 0.220, 0.400, 0.250),
#         "tmax" = list(0.130, 0.180, 0.300, 0.600, 0.300),
#         "roi" = list(
#           c("O1", "O2", "Oz", "PO3", "PO4", "PO7", "PO8"),
#           c("O1", "O2", "Oz", "PO3", "PO4", "PO7", "PO8"),
#           c("Fz", "Cz"),
#           c("Fz", "Cz", "CPz", "Pz"),
#           c("CP3", "CP4", "P5", "P3", "Pz", "P4", "P6", "PO3", "PO4"))
#         ),
#    average_by = c("test_part/erp2_out/naming",
#                   "test_part/erp2_out/frequency",
#                   "test_part/erp2_out/naming/frequency"),
#    n_jobs=3
#    )


# Since the RIDE implementation into the pipeline is not launched yet, 
# we copmpute it directly in python

#import os
#import sys

from pipeline import group_pipeline
#sys.path.insert(1, os.path.join(sys.path[0], '..'))
from pyprojroot import here

# all participants - memory task
res=group_pipeline(
    vhdr_files=here("data/eeg_VP1/raw/memory/"),
    log_files=here("data/raw/log/with_practice_for_EEG_memory/"),
    output_dir=here("data/eeg_VP1/export/memory_ERP2_ride"),
   #report_dir=here("data/eeg_VP1/export/memory_ERP2_ride/qc_reports"),
   #epochs_dir=
   #clean_dir=
   besa_files=None,
   ica_method="fastica",
   ica_n_components=0.99,
   highpass_freq=0.1,
   lowpass_freq=40.0,
   reject_peak_to_peak=200.0,
    triggers=[221, 222], # 222
    triggers_column = "bedingung_trigger_erp2",
    epochs_tmin = -0.2,
    epochs_tmax = 3.2,
    baseline = [-0.2, 0.0],
    bad_channels = { "VP0301" : ["T7"],
                     "VP0330" : ["TP10"]},
    perform_ride=True,
    ride_condition_column='bedingung_trigger_erp2', # Kombination aus 2 Spalten?
    ride_rt_column='vkRT',
    ride_s_twd=(0.0, 0.6),
    ride_r_twd=(-0.3, 0.3),
    ride_epochs_tmin_after_ride=-0.2,
    ride_epochs_tmax_after_ride=0.8,
    ride_reject_peak_to_peak= 200.0,
   # skip_log_conditions = {'bedingung_trigger_erp2': [221]},
    components = { 'name' : ['P1','P1_local', 'P1_local_c', 'N1', 'N1_local', 'N2', 'N2_local', 'P300', 'posP2', 'posP2_local', 'posP2_local_c', 'EPN_ex', 'LPP_ex'],
                    'tmin' : [0.080, 0.091, 0.093, 0.130, 0.143, 0.220, 0.204, 0.400, 0.250, 0.219, 0.209, 0.200, 0.400],
                    'tmax' : [0.130, 0.141, 0.143, 0.180, 0.193, 0.300, 0.284, 0.600, 0.300, 0.269, 0.259, 0.300, 0.600],
                    'roi' : [['O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8'],
                             ['O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8'],
                             ['O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8'],
                            ['O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8'],
                            ['O1', 'O2', 'Oz', 'PO3', 'PO4', 'PO7', 'PO8'],
                            ['Fz', 'Cz'],
                            ['Fz', 'Cz'],
                            ['Fz', 'Cz', 'CPz', 'Pz'],
                            ['CP3', 'CP4', 'P5', 'P3', 'Pz', 'P4', 'P6', 'PO3', 'PO4'],
                            ['CP3', 'CP4', 'P5', 'P3', 'Pz', 'P4', 'P6', 'PO3', 'PO4'],
                            ['CP3', 'CP4', 'P5', 'P3', 'Pz', 'P4', 'P6', 'PO3', 'PO4'],
                            ['PO7', 'PO8', 'PO9', 'PO10', 'TP9', 'TP10'],
                            ['Pz', 'CPz', 'POz', 'P3', 'P4']]},
    average_by = { #"test_trials": "test_part=='practice-memory'",
                #"nan_trials": "vkRT.isna()",
                #"0_trials": "vkRT==0",
                "simple_naming": "test_part=='memory-task' & bedingung_trigger_erp2==221 & erp2_out=='in'",
                "memory_naming": "test_part=='memory-task' & bedingung_trigger_erp2==222 & erp2_out=='in'",
                "low_frequency": "test_part=='memory-task' & frequency=='low' & erp2_out=='in'",
                "high_frequency": "test_part=='memory-task' & frequency=='high' & erp2_out=='in'",
                "simple_naming_low": "test_part=='memory-task' & bedingung_trigger_erp2==221 & VLcodierung==250 & frequency=='low' & erp2_out=='in'",
                "simple_naming_high": "test_part=='memory-task' & bedingung_trigger_erp2==221 & VLcodierung==250 & frequency=='high' & erp2_out=='in'",
                "memory_naming_low": "test_part=='memory-task' & bedingung_trigger_erp2==222 & VLcodierung==251 & frequency=='low' & erp2_out=='in'",
                "memory_naming_high": "test_part=='memory-task' & bedingung_trigger_erp2==222 & VLcodierung==251 & frequency=='high' & erp2_out=='in'",
               },
    n_jobs=2
        # add rois and frequencies
            # ["test_part/erp2_out/naming"] #,
                 #  "test_part/erp2_out/frequency",
                 #  "test_part/erp2_out/naming/frequency"]
    ) 
