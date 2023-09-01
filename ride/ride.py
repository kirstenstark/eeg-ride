import numpy as np
from pyprojroot.here import here
from scipy.fft import fft, ifft
from scipy.interpolate import CubicSpline
from scipy.io import loadmat, savemat

# import matplotlib.pyplot as plt # for plotting

# from tests import generate_sine_data, plot_input_data

# # Generate simulated sine data
# data, rt = generate_sine_data()
# plot_input_data(data)

# mat_file = here('data/sine_data.mat')
# mat_dict = {'data': data}
# savemat(mat_file, mat_dict)

# load exemplary EEG data

mat_file = here('matlab/matlab_example/mat/Vp0001_rep1_close.mat')
data_dict = loadmat(str(mat_file), mat_dtype=True)
data = data_dict['data']
data = data.astype('float64')
rt = data_dict['rt']
rt = rt.astype('float64')

# define cfg file
cfg = {'samp_interval': 2,
       'epoch_twd': np.array([-100, 1198]),
       'comp': {'name': ['s', 'r'],
                'twd': [[0, 600], [-300, 300]],
                'latency': [0, rt]}}
cfg['comp_num'] = len(cfg['comp']['name'])
cfg['rwd'] = 200
cfg['re_samp'] = cfg['samp_interval']
cfg['bd'] = 0.2  # Alpha value for Tukey window
cfg['prg'] = 1  # Print progress to console
cfg['bl'] = 200  # Baseline length

# start RIDE correction
assert len(cfg['comp']['name']) > 1, 'At least two components are required'

# TODO: make sure order is always s, (c), r

cfg0 = cfg.copy()

# section 1
d1, d2, d3 = data.shape
epoch_length = d1
erp = data.mean(axis=2)
results = {'latency0': cfg['comp']['latency']}


def round_like_matlab(x):
    """Round to nearest integer, like MATLAB's round() function."""

    return np.trunc(x + np.copysign(0.5, x)).astype(int)


rs = cfg['re_samp'] / cfg['samp_interval']
data = data[round_like_matlab(np.linspace(0, d1 - 1, int(d1 / rs))), :, :]
d1, d2, d3 = data.shape  # New size after down samping


for j in range(cfg['comp_num']):

    if isinstance(cfg['comp']['latency'][j], int):
        int_value = cfg['comp']['latency'][j]
        print(f'WARNING: Extending integer latency {int_value} to a vector of {int_value}s (one per trial)')
        cfg['comp']['latency'][j] = np.array([[int_value]] * d3)

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
c_sl = c_l.copy()

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


def ride_detrend(data, twd):

    # # For debugging
    # data = temp.copy()
    # twd = np.array([51, 110, 288, 348])

    index = np.isnan(data)
    data[index] = 0.0
    d1, d2 = data.shape
    d3 = 1  # MATLAB uses `d1, d2, d3 = data.shape` even though `data` is 2D
    d = [(twd[3] + twd[2]) / 2 - (twd[1] + twd[0]) / 2]

    temp0 = np.tile(np.arange(d1) + 1, [d2, d3]).T

    a = data[np.arange(twd[0], twd[1] + 1), :].mean(axis=0, keepdims=True)
    b = data[np.arange(twd[2], twd[3] + 1), :].mean(axis=0, keepdims=True)
    temp = (b - a) / d
    data = data - temp0 * temp[np.zeros(d1, dtype=int), :]

    temp = data[np.arange(twd[0], twd[1] + 1), :].mean(axis=0, keepdims=True)
    data = data - temp[np.zeros(d1, dtype=int), :]

    f = data.copy()
    f[index] = np.nan

    return f


def median_2d(x):
    # !! Median is probably biased if nan trials are present b/c they are set to zero
    # TODO: simulate this and see if it's a problem
    x = x.copy()
    x[np.isnan(x)] = 0.0  # zeros padding
    s = x.shape
    x = np.sort(x, axis=0)
    y = x.flatten(order='F')[(round_like_matlab(
        s[0]/2) + np.arange(0, s[1])*s[0])-1]
    return y


def RIDE_tukey(n, r):
    n = int(n)
    t = np.linspace(0, 1, n)
    # Define period of the taper as 1/2 period of a sine wave.
    per = r/2
    tl = int(np.floor(per*(n-1))+1)
    th = n-tl+1
    # Window is defined in three sections: taper, constant, taper
    f = np.concatenate([(1+np.cos(np.pi/per*(t[np.arange(tl)] - per)))/2,
                       np.ones((th-tl-1)),
                       (1+np.cos(np.pi/per*(t[th-1:] - 1 + per)))/2])

    if n == 1:
        f = 1

    return f

def mean_nan(data, dim):
    data = np.nan_to_num(data)
    f = np.mean(data, dim)
    return f

def ride_iter(data, cfg):

    # # For debugging
    # cfg = cfg1
    # data = data[:, 61, :]

    d1, d2 = data.shape
    bd = cfg['bd']

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
    for c in range(cfg['comp_num']):
        if cfg['comp']['name'][c] == 'r':
            trend_i = c-1
    if trend_i == cfg['comp_num']:
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i)])
    if trend_i == cfg['comp_num']-2 : 
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i), np.array([cfg['comp_num']-1])])
    if trend_i == 0 : 
        stream_flow = np.array([1,0])

    l1 = np.zeros((cfg['inner_iter'], cfg['comp_num']))
    stop = 0
    stop_c = np.zeros((cfg['comp_num'], 1))
    com_old = com_c.copy()

    # for iter in np.arange(1, 2): # For testing purposes, only runs iter = 1
    for iter in np.arange(cfg['inner_iter']):

        if iter + 1 == cfg['inner_iter']:
            stop = 1

        # track the convergence

        # !!!! CONTINUE HERE, AFTER FINISHING FIRST ITERATION OF FOR LOOP
        # TO DO: add if-loop for iter > 1

        if iter > 1: 
            for c in np.arange(cfg['comp_num']):
                l1[iter-2,c]=np.sum(np.abs(com_c[:,c] - com_old[:,c]))
                if iter > 2:
                    if l1[iter-3,c] - l1[iter-2,c] < 0.01*(l1[0,c] - l1[1,c]): # convergence is now defined as 0.01
                        stop_c[c] = 1
                    if l1[iter-3,c] == l1[iter-2,c]:
                        stop_c[c] = 1
            if cfg['comp_num'] == 1:
                stop_c[c] = 1
        for c in np.arange(cfg['comp_num']):
            com_old[:,c] = com_c[:,c]
            # decision of the termination of the interation of each component
            

        
        #for c in np.arange(cfg['comp_num']):
            #com_old[:,c]=com_c[:,c] # decision of the termination of iteration of each component
        ## CAUTION: for-loop seems unnecessary

        # com_old = com_c.copy()

        for c in stream_flow:
            if stop_c[c] == 0:
                temp = data.copy()
                for j in np.arange(cfg['comp_num']):
                    if j != c:
                        temp = temp-com_c1[:, :, j]
                residue = temp.copy()
                temp = np.empty((length_c[c], d2))
                temp[:] = np.nan
                for j in np.arange(d2):
                    temp[np.arange(-cfg['comp']['latency'][c][j]+max_latency[c],d1-cfg['comp']['latency'][c][j]+max_latency[c]),j] = residue[:,j]

                temp = ride_detrend(
                    temp,
                    np.array([0,
                              np.fix((cfg['comp']['twd'][c][1] - cfg['comp']['twd'][c][0]) * bd),
                              np.fix((cfg['comp']['twd'][c][1] - cfg['comp']['twd'][c][0]) * (1 - bd)) - 1,
                              cfg['comp']['twd'][c][1] - cfg['comp']['twd'][c][0] - 1], dtype=int) \
                                + max_latency[c] + int(cfg['comp']['twd'][c][0]))
                
                temp0 = median_2d(temp.T)
                temp0[np.concatenate([np.arange(cfg['comp']['twd'][c][0], dtype=int) + max_latency[c],
                                      np.arange(cfg['comp']['twd'][c][1] + max_latency[c]-1, length_c[c], dtype=int)])]=0
                temp0[np.arange(cfg['comp']['twd'][c][0]+ max_latency[c]-1, cfg['comp']['twd'][c][1] + max_latency[c], dtype=int)] = \
                    temp0[np.arange(cfg['comp']['twd'][c][0]+ max_latency[c]-1, cfg['comp']['twd'][c][1] + max_latency[c], 
                                    dtype=int)] * RIDE_tukey(cfg['comp']['twd'][c][1] - cfg['comp']['twd'][c][0] + 1, bd*2)
     
                temp1 = np.repeat(temp0[:,np.newaxis],d2, axis=1)
                temp1[np.isnan(temp)] = np.nan

                # TODO: !!!! IMPORTANT: Next line is only for testing purposes, delete later!!
                # temp1[:,44] = np.nan     # # TODO: !!!! IMPORTANT: This line is only for testing purposes, delete later!!
                temp1 = np.reshape(temp1.T[~np.isnan(temp1.T)], (d1, d2), order='F')
                com_c[:,c] = temp0[np.arange(max_latency[c], max_latency[c]+d1)]
                com_c1[:,:,c] = temp1
            # End of if-loop
        # End of stream-flow for loop

        if stop_c.mean() == 1:
            stop = 1
        if stop == 1:
            for c in stream_flow:

                temp = data.copy()
                for j in np.arange(cfg['comp_num']):
                    if j != c:
                        temp = temp - com_c1[:,:,j]

                residue = temp.copy()
                temp = np.empty((length_c[c], d2))
                temp[:] = np.nan

                for j in np.arange(d2):
                    temp[np.arange(-cfg['comp']['latency'][c][j]+max_latency[c],
                                   d1-cfg['comp']['latency'][c][j]+max_latency[c]),
                         j] = residue[:,j]

                if cfg['final'] == 1:
                    temp[np.isnan(temp)] = 0
                    amp_c[:,c] = np.mean(
                        np.repeat(com_c[np.arange(cfg["comp"]["twd"][c][0] - 1, cfg["comp"]["twd"][c][1], dtype=int), c, np.newaxis], d2, axis=1) * temp[max_latency[c] + np.arange(cfg["comp"]["twd"][c][0] - 1, cfg["comp"]["twd"][c][1], dtype=int), :], axis=0)

            break
        # end of inner iter loop
    
    # for last iteration
    if cfg['final'] == 1:
        # release time window function and detrending
        # allocate trend to the last C component
        for c in stream_flow:
            temp = data0.copy()
            for j in np.arange(cfg['comp_num']):
                if j != c:
                    temp = temp - com_c1[:,:,j]
                    # TODO: Probably only for microsaccades, temp-com_ms1 if 'ms' and 'var' exist
            residue = temp.copy()
            temp = np.empty((length_c[c], d2))
            temp[:] = np.nan

            for j in np.arange(d2):
                temp[np.arange(-cfg['comp']['latency'][c][j]+max_latency[c],
                               d1-cfg['comp']['latency'][c][j]+max_latency[c]),
                               j] = residue[:,j]
            temp0 = mean_nan(temp,1)

            if c==stream_flow[0]:
                temp0[np.arange(cfg['comp']['twd'][c][0] + max_latency[c],dtype=int)]=0.0
                tem = RIDE_tukey(cfg['comp']['twd'][c][1] - cfg['comp']['twd'][c][0] + 1, bd*2)
                tem = tem[np.arange(np.fix(len(tem)/2).astype(int))]
                temp0[np.arange(cfg['comp']['twd'][c][0] + max_latency[c]-1, 
                                cfg['comp']['twd'][c][0] + max_latency[c] + len(tem)-1, dtype=int)] = \
                                    temp0[np.arange(cfg['comp']['twd'][c][0] + max_latency[c]-1,\
                                                    cfg['comp']['twd'][c][0] + max_latency[c] + len(tem)-1, dtype=int)] * tem
                
            temp1 = np.repeat(temp0[:,np.newaxis],d2, axis=1)
            temp1[np.isnan(temp)] = np.nan
            temp1 = np.reshape(temp1.T[~np.isnan(temp1.T)], (d1, d2), order='F')
            com_c[:,c] = temp0[np.arange(max_latency[c], max_latency[c]+d1)]
            com_c1[:,:,c] = temp1.copy()
        
        c = stream_flow[0]
        temp = data0.copy()
        for j in np.arange(cfg['comp_num']):
            if j != c:
                temp = temp - com_c1[:,:,j]
                # TODO: Probably only for microsaccades, temp-com_ms1 if 'ms' and 'var' exist
        residue=temp.copy()
        temp = np.empty((length_c[c], d2))
        temp[:] = np.nan
        for j in np.arange(d2):
            temp[np.arange(-cfg['comp']['latency'][c][j]+max_latency[c],
                           d1-cfg['comp']['latency'][c][j]+max_latency[c]),
                           j] = residue[:,j]
        temp0 = mean_nan(temp,1)
        temp1 = np.repeat(temp0[:,np.newaxis],d2, axis=1)
        temp1[np.isnan(temp)] = np.nan
        temp1 = np.reshape(temp1.T[~np.isnan(temp1.T)], (d1, d2), order='F')
        com_c[:,c] = temp0[np.arange(max_latency[c], max_latency[c]+d1)]
        com_c1[:,:,c] = temp1.copy()

    # end of last iteration
    results={}
    results['amp'] = amp_c
    results['comp'] = com_c
    results['comp1'] = np.mean(com_c1,1)
    results['iter'] = iter
    results['l1'] = l1
    results['trend_c'] = stream_flow[0]

    return results

def interp2d(data, x1, x2):

    temp = np.empty((len(x2), data.shape[1]))

    for j in np.arange(data.shape[1]):
        cs = CubicSpline(x1, data[:, j], axis=0)
        temp[:, j] = cs(x2)
    
    return temp

def baseline(x, a=None):

    x[np.isnan(x)] = 0.0

    if a is None:
        temp = np.mean(x, axis=0)
        f = x - temp[np.newaxis, :]
    else:
        temp = np.mean(x[a, :], axis=0)
        f = x - temp[np.newaxis, :]

    return f


amp = np.zeros((d3, d2, cfg['comp_num']))
for c in range(d2):

    if cfg['prg'] == 1:
        print(f'Processing electrode #{c}')

    rst = ride_iter(np.squeeze(data[:, c, :]), cfg1)

    c_l[:, :, c] = rst['comp']
    c_sl[:, :, c] = rst['comp1']

    if stop == 1:
        amp[:, c, :] = rst['amp']
    
    if cfg['prg'] == 1:
        print(f'Took {rst["iter"] + 1} iterations')

comp = c_l.transpose((0, 2, 1))
comp1 = c_sl.transpose((0, 2, 1))

results['erp_new'] = 0
results['residue'] = erp

bl_wd = np.arange(-cfg['epoch_twd'][0]/cfg['samp_interval'],
                  -cfg['epoch_twd'][0]/cfg['samp_interval']+cfg['bl']/cfg['samp_interval'],
                  dtype=int)


component = np.zeros((d1, d2, cfg['comp_num']))
component1 = np.zeros((d1, d2, cfg['comp_num']))

for j in np.arange(cfg['comp_num']):
    # The MATLAB version explicitly requests "spline" as the interpolation method
    # Our Python function `interp2d` only performs this spline interpolation,
    # with no support for any other method (unlike the MATLAB function)
    component[:, : , j] = interp2d(comp[:, :, j],
                                   np.round(np.linspace(0, epoch_length, d1, endpoint=False)),
                                   np.arange(0, epoch_length))
    component[:, :, j] = baseline(component[:, :, j],bl_wd)
    component1[:, :, j] = interp2d(comp1[:, :, j],
                                    np.round(np.linspace(0, epoch_length, d1, endpoint=False)),
                                    np.arange(0, epoch_length))
    component1[:, :, j] = baseline(component1[:, :, j],bl_wd)
    results['residue'] = results['residue'] - component1[:, :, j]
    results[cfg['comp']['name'][j]] = component[:, :, j]
    results[cfg['comp']['name'][j]+'_sl'] = component1[:, :, j]
    results['latency_'+cfg['comp']['name'][j]] = cfg['comp']['latency'][j]*cfg['re_samp']
    results['amp_'+cfg['comp']['name'][j]] = amp[:, :, j]

results[cfg['comp']['name'][0]] = baseline(results[cfg['comp']['name'][0]],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)
results[cfg['comp']['name'][0]+'_sl']= baseline(results[cfg['comp']['name'][0]+'_sl'],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)

results[cfg['comp']['name'][rst['trend_c']]] = baseline(results[cfg['comp']['name'][rst['trend_c']]] + results['residue'],bl_wd)
results[cfg['comp']['name'][rst['trend_c']]+'_sl'] = baseline(results[cfg['comp']['name'][rst['trend_c']]+'_sl'] + results['residue'],bl_wd)

if  cfg['comp_num'] == 1:
    bl_wd = np.arange(-cfg['epoch_twd'][0]/cfg['samp_interval'], dtype='int')
    results[cfg['comp']['name'][0]] = baseline(results[cfg['comp']['name'][0]] + results['residue'],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)
    results[cfg['comp']['name'][0]+'_sl'] = baseline(results[cfg['comp']['name'][0]+'_sl'] + results['residue'],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)

for j in np.arange(cfg['comp_num']):
    results['erp_new'] = results['erp_new'] + results[cfg['comp']['name'][j]]

results['cfg'] = cfg0.copy()
