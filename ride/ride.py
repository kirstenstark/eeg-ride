import numpy as np

from .baseline import baseline
from .helpers import round_like_matlab
from .interp2d import interp2d
from .iter import ride_iter


def ride_call(data, cfg):

    # # Load example data
    # from pyprojroot.here import here
    # mat_file = here('matlab/matlab_example/mat/Vp0008_rep1_distant.mat')
    # data_dict = loadmat(str(mat_file), mat_dtype=True)
    # data = data_dict['data']
    # data = data.astype('float64')
    # rt = data_dict['rt']
    # rt = rt.astype('float64')

    # start RIDE correction
    assert len(cfg['comp']['name']) > 1, 'At least two components are required'

    # TODO: make sure order is always s, (c), r

    cfg0 = cfg.copy()

    # section 1
    d1, d2, d3 = data.shape
    epoch_length = d1
    erp = data.mean(axis=2)
    results = {'erp': erp,
            'latency0': cfg['comp']['latency'].copy()}

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

    return results
