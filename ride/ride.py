from pyprojroot.here import here
from scipy.io import savemat, loadmat
from tests import generate_sine_data, plot_input_data
import numpy as np
 
# # Generate simulated sine data
# data, rt = generate_sine_data()
# plot_input_data(data)

# mat_file = here('data/sine_data.mat')
# mat_dict = {'data': data}
# savemat(mat_file, mat_dict)

# load exemplary EEG data

mat_file = here('matlab/matlab_example/mat/Vp0001_rep1_close.mat')
data_dict = loadmat(str(mat_file))
data = data_dict['data']
rt = data_dict['rt']

# define cfg file
cfg = {'samp_interval': 2,
       'epoch_twd': np.array([-100, 1198]),
       'comp': {'name': ['s', 'r'],
                 'twd': [[0, 600], [-300, 300]],
                 'latency': [0, rt]}}
cfg['comp_num'] = len(cfg['comp']['name'])
cfg['rwd'] = 200
cfg['re_samp'] = cfg['samp_interval']

# start RIDE correction
assert len(cfg['comp']['name']) > 1, 'At least two components are required'

cfg0 = cfg.copy()

## section 1
d1, d2, d3 = data.shape
epoch_length = d1
erp = data.mean(axis=2)
results={'latency0': cfg['comp']['latency']}

# TODO: Do downsampling here and re-obtain d1, d2, d3


def round_like_matlab(x):
    """Round to nearest integer, like MATLAB's round() function."""

    return np.trunc(x + np.copysign(0.5, x))


for j in range(cfg['comp_num']):
    
    if isinstance(cfg['comp']['latency'][j], int):
        int_value = cfg['comp']['latency'][j]
        print(f'WARNING: Extending integer latency {int_value} to a vector ' +
              f'of {int_value}s (one per trial)')
        cfg['comp']['latency'][j] = np.array([int_value] * d3)
    
    if cfg['comp']['name'][j] == 'r':
        cfg['comp']['twd'][j] = cfg['comp']['twd'][j] + np.median(results['latency0'][j])
        cfg['comp']['twd'][j][cfg['comp']['twd'][j] < cfg['rwd']] = cfg['rwd']
        cfg['comp']['twd'][j][cfg['comp']['twd'][j] > cfg['epoch_twd'][1]] = cfg['epoch_twd'][1]

    cfg['comp']['latency'][j] = cfg['comp']['latency'][j] / cfg['re_samp']
    cfg['comp']['latency'][j] = round_like_matlab(cfg['comp']['latency'][j]-np.median(cfg['comp']['latency'][j])) # Still floats, might need to be ints
    
    cfg['comp']['twd'][j] = np.fix((cfg['comp']['twd'][j] - cfg['epoch_twd'][0])/cfg['re_samp'])+[1,-1]

