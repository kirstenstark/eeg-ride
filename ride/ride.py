import numpy as np
from pyprojroot.here import here
from scipy.fft import fft, ifft
from scipy.io import loadmat, savemat
# from tests import generate_sine_data, plot_input_data

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

# TODO: make sure order is always s, (c), r

cfg0 = cfg.copy()

# section 1
d1, d2, d3 = data.shape
epoch_length = d1
erp = data.mean(axis=2)
results = {'latency0': cfg['comp']['latency']}

# TODO: Do downsampling here and re-obtain d1, d2, d3


def round_like_matlab(x):
    """Round to nearest integer, like MATLAB's round() function."""

    return np.trunc(x + np.copysign(0.5, x)).astype(int)


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
    cfg['comp']['latency'][j] = round_like_matlab(cfg['comp']['latency'][j]-np.median(cfg['comp']['latency'][j])) 

    cfg['comp']['twd'][j] = np.fix((cfg['comp']['twd'][j] - cfg['epoch_twd'][0])/cfg['re_samp'])+[1, -1]


stop = 1

# TODO: Add "outer iteration" loop around here if there are one or more C components

cfg1 = cfg.copy()
cfg1['final'] = stop
cfg1['inner_iter'] = 100

c_l = np.zeros((d1, cfg['comp_num'], d2))
c_sl = c_l

# for c in range(d2):
#     rst = ride_iter(np.squeeze(data[:, c, :]), cfg1)


def filtering20(x, a, b):

    f = x
    for j in range(x.shape[1]):
        f[:, j] = filtering10(x[:, j], a, b)

    return f


def filtering10(x, a, b):

    # The MATLAB code does x = x(:) here but this doesn't seem to change anything
    x0 = x.mean()
    x = x - x0
    n = 10 * len(x)
    Y = fft(x, n)

    b = int(b)
    H = np.concatenate([np.zeros((1, a)),
                        np.expand_dims(Y[a: b + 1], axis=0),
                        np.zeros((1, n - b * 2 - 1)),
                        np.expand_dims(Y[n - b: n - a + 1], axis=0),
                        np.zeros((1, a - 1))], axis=1)
    H = H.transpose()

    f = ifft(H, axis=0).real
    f = np.squeeze(f)

    return f[:len(x)]


def ride_iter(data, cfg):
    
    # These are the inputs that get passed to the function
    # Remove these two lines once the function is complete!
    cfg = cfg1
    data = data[:, 61, :]

    d1, d2 = data.shape

    data0 = data.copy()

    data = filtering20(data, 1, np.fix(10 * 20 * cfg['re_samp'] * d1 / 1000))
    # TODO: Compare with MNE filter

    max_latency = np.zeros((cfg['comp_num']), dtype=int)
    min_latency = np.zeros((cfg['comp_num']), dtype=int)
    length_c = np.zeros((cfg['comp_num']), dtype=int)
    for j in range(cfg['comp_num']):
        max_latency[j] = cfg['comp']['latency'][j].max()
        min_latency[j] = cfg['comp']['latency'][j].min()
        length_c[j] = d1 + max_latency[j] - min_latency[j]

    com_c = np.zeros((d1, cfg['comp_num']))
    com_c1 = np.zeros((d1, d2, cfg['comp_num']))
    amp_c = np.zeros((d2, cfg['comp_num']))

    trend_i = cfg['comp_num']-1
    stream_flow = np.arange(cfg['comp_num'])
    for c in range(cfg['comp_num']) :
        if cfg['comp']['name'][c] == 'r':
            trend_i = c-1
    if trend_i == cfg['comp_num'] : 
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i)])
    if trend_i == cfg['comp_num']-2 : 
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i), np.array([cfg['comp_num']-1])])
    if trend_i == 0 : 
        stream_flow = np.array([1,0])

    l1 = np.zeros((cfg['inner_iter'],cfg['comp_num']))
    stop = 0
    stop_c = np.zeros((cfg['comp_num'],1))
    com_old=com_c.copy()

    for iter in np.arange(cfg['inner_iter']) :
        
        if iter + 1 == cfg['inner_iter']:
            stop = 1
        
        # track the convergence

        # !!!! CONTINUE HERE, AFTER FINISHING FIRST ITERATION OF FOR LOOP
        # TO DO: add if-loop for iter > 1

        #if iter > 1: 
          #  for c in np.arange(cfg['comp_num']):
           #     l1[iter-2,c]=np.sum()
                # continue here
        
        #for c in np.arange(cfg['comp_num']):
            #com_old[:,c]=com_c[:,c] # decision of the termination of iteration of each component
        ## CAUTION: for-loop seems unnecessary
        com_old=com_c.copy()

        for c in stream_flow: 
            if stop_c[c] == 0:
                temp = data.copy()
                for j in np.arange(cfg['comp_num']):
                    if j != c:
                        temp = temp-com_c1[:,:,j]
                residue = temp.copy()
                temp = np.empty((length_c[c],d2))
                temp[:] = np.nan
                for j in np.arange(d2):
                    temp[np.arange(-cfg['comp']['latency'][c][j]+max_latency[c],d1-cfg['comp']['latency'][c][j]+max_latency[c]),j] = residue[:,j]

                 ### CONTINUE HERE ###
