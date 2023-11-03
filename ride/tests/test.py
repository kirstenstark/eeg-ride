import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from pyprojroot.here import here
from scipy.io import loadmat

sys.path.insert(1, os.path.join(sys.path[0], '..'))

from ride.cfg import RideCfg
from ride.ride import ride_call

mat_file = here('matlab/matlab_example/mat/Vp0001_rep1_close.mat')
data_dict = loadmat(str(mat_file), mat_dtype=True)
data = data_dict['data']
data = data.astype('float64')
rt = data_dict['rt']
rt = rt.astype('float64')


cfg = RideCfg(comp_name=['s', 'r'],
              comp_twd=[[0, 600], [-300, 300]],
              comp_latency=[0, rt],
              sfreq=500)

results = ride_call(data, cfg)

n_comps = len(results.comps)
ncols = 1 + n_comps
fig, axs = plt.subplots(1, ncols, figsize=(10, 3))
x = np.linspace(cfg.epoch_twd[0], cfg.epoch_twd[1], data.shape[0])

axs[0].plot(x, results.erp, linewidth=0.7)
axs[0].set_title('erp')
axs[0].set_ylim([-20.0, 20.0])
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('Amplitude (µV)')

for ax, key in zip(axs[1:], results.comps):
    ax.plot(x, results.comps[key], linewidth=0.7)
    ax.set_title(key)
    ax.set_ylim([-20.0, 20.0])
    ax.set_xlabel('Time (ms)')

fig.show()
