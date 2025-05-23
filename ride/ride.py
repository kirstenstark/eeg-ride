from warnings import warn

import numpy as np
from mne import Epochs, EpochsArray

from .baseline import baseline
from .helpers import round_like_matlab
from .interp2d import interp2d
from .iter import ride_iter
from .results import RideResults

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
    assert len(cfg.comp_name) > 1, 'At least two components are required'

    # TODO: make sure order is always s, (c), r

    cfg0 = cfg.copy()

    # If necessary, extract Numpy arra from MNE epochs object
    if isinstance(data, (Epochs, EpochsArray)):
        data = data.get_data(picks='eeg')
        data = np.swapaxes(data, 0, 2)

    # section 1
    d1, d2, d3 = data.shape
    epoch_length = d1
    erp = data.mean(axis=2)
    results = RideResults(erp=erp, latency0=cfg.comp_latency.copy(), latency_all=cfg.comp_latency.copy())

    rs = cfg.re_samp / cfg.samp_interval
    data = data[round_like_matlab(np.linspace(0, d1 - 1, int(d1 / rs))), :, :]
    d1, d2, d3 = data.shape  # New size after down samping
    
    nan_ixs = set()

    # define nan_ixs (and RT = 0 trials)
    for j in range(cfg.comp_num):

        if isinstance(cfg.comp_latency[j], (int, float)):
            int_value = cfg.comp_latency[j]
            if cfg.prg == 1:
                print(f'WARNING: Extending integer latency {int_value} to a vector of {int_value}s (one per trial)')
            cfg.comp_latency[j] = np.array([[int_value]] * d3)

        if not isinstance(cfg.comp_latency[j], str):
            cfg.comp_latency[j] = np.array(cfg.comp_latency[j])
            
            nan_ixs_comp = np.where(np.isnan(cfg.comp_latency[j]))[0]
            nan_ixs.update(nan_ixs_comp)

            if cfg.comp_name[j] == 'r':
                nan_ixs_comp = np.where(cfg.comp_latency[j]==0)[0]
                nan_ixs.update(nan_ixs_comp)
    
    nan_ixs = list(nan_ixs)
    if len(nan_ixs) > 0:
        warn(f'Trials {nan_ixs} will not be used for RIDE estimation because their latency is NaN or zero.')
        data = np.delete(data, nan_ixs, axis=2)
        d1, d2, d3 = data.shape  # New size after nan removal
        for j in range(cfg.comp_num):
            cfg.comp_latency[j] = np.delete(cfg.comp_latency[j], nan_ixs, axis=0)
            if cfg.comp_name[j] == 'r':
                results.latency0[j] = np.delete(results.latency0[j], nan_ixs, axis=0)

    # Resample data
    for j in range(cfg.comp_num):
        if cfg.comp_name[j] == 'r':
            cfg.comp_twd_samp[j] = cfg.comp_twd_samp[j] + np.median(results.latency0[j])
            cfg.comp_twd_samp[j][cfg.comp_twd_samp[j] < cfg.rwd] = cfg.rwd
            cfg.comp_twd_samp[j][cfg.comp_twd_samp[j] > cfg.epoch_twd[1]] = cfg.epoch_twd[1]

        if not isinstance(cfg.comp_latency[j], str):
            cfg.comp_latency[j] = cfg.comp_latency[j] / cfg.re_samp
            cfg.comp_latency[j] = round_like_matlab(cfg.comp_latency[j]-np.median(cfg.comp_latency[j]))

        cfg.comp_twd_samp[j] = np.array(np.fix((cfg.comp_twd_samp[j] - cfg.epoch_twd[0])/cfg.re_samp)+[1, -1], dtype=int)

    for j in range(cfg.comp_num):
        if cfg.dur[j] is not None:
            cfg.dur[j] = int(np.fix(cfg.dur[j] / cfg.re_samp))
        else:
            cfg.dur[j] = round((cfg.comp_twd_samp[j][1] - cfg.comp_twd_samp[j][0]) / 2)

    # Initial estimation of the latency of C component
# %     for section = 1:1%initial estimation of the latency of C---------------------------------------
# %     for initial_c = 1:1 
# %             n_of_c = 0;c_i = 0;
# %             for j = 1:cfg.comp_num
# %                 if ischar(cfg.comp.latency{j})
# %                         if cfg.prg == 1 disp(['woody_for_',cfg.comp.name{j}]);end
# %                         n_of_c = n_of_c + 1;c_i(n_of_c) = j;
# %                         temp = 1:d2;temp1 = cfg.comp.twd{j};
# %                             if isfield(cfg,'template')
# %                                 if strcmp(cfg.template.method,'g_mean')
# %                                     cfg.temp = template(temp1(1):temp1(2),temp);
# %                                 end
# %                                 if isfield(cfg.template,'chan')
# %                                     temp = cfg.template.chan;
# %                                     if isfield(cfg.template,'hann_amp')
# %                                         cfg.template.hann_amp = cfg.template.hann_amp(cfg.template.chan);
# %                                     end
# %                                 end
# %                             end
# %                             %-------using Woody's method by default
# %                         cfg.comp.latency{j} = woody(data(temp1(1):temp1(2),temp,:),cfg,cfg.dur{j});%
                        
# %                 end
# %             end
# %     end
# % end
# Translate the entire section to Python
    n_of_c = 0
    c_i = []
    for j in range(cfg.comp_num):
        if isinstance(cfg.comp_latency[j], str): # ToDo: Change latency of C component from 'unknown' to xx
            if cfg.prg == 1:
                print(f'woody_for_{cfg.comp_name[j]}')
            n_of_c += 1
            c_i.append(j)
            temp = np.arange(d2)
            temp1 = cfg.comp_twd_samp[j]
            # ToDo: Check code once 'template' is implemented
            if hasattr(cfg, 'template'):
                if cfg.template.method == 'g_mean':
                    cfg.temp = cfg.template[temp1[0]:temp1[1], temp]
                if hasattr(cfg.template, 'chan'):
                    temp = cfg.template.chan
                    if hasattr(cfg.template, 'hann_amp'):
                        cfg.template.hann_amp = cfg.template.hann_amp[cfg.template.chan]
            cfg.comp_latency[j] = woody(data[(temp1[0]-1):temp1[1], temp, :],
                                        cfg, cfg.dur[j])

    stop = 1

    latency_i = np.zeros(n_of_c)
    for j in range(n_of_c): # Track latency evolution of C components
        # TODO: Continue translating outer loop here
        # latency_i[j] = cfg.comp.latency[c_i[j]]
        # latency_i[j] = latency_i[j](:);
    # %     l_change(:,j) = ones(d3,1);%track for evolution of the latency in order to terminate the iteration
    # %     c_change(:,j) = ones(d3,1);%track for evolution of the correlation in order to terminate the iteration
    # % end

    # % if cfg.prg == 1 disp('RIDE decomposition: ');end
    # outer_iter = 4;if n_of_c == 0 outer_iter = 1;end

    cfg1 = cfg.copy()
    cfg1.final = stop
    cfg1.inner_iter = 100

    c_l = np.zeros((d1, cfg.comp_num, d2))
    c_sl = c_l.copy()

    amp = np.zeros((d3, d2, cfg.comp_num))
    for c in range(d2):

        if cfg.prg == 1:
            print(f'Processing electrode #{c + 1}')

        rst = ride_iter(np.squeeze(data[:, c, :]), cfg1)

        c_l[:, :, c] = rst['comp']
        c_sl[:, :, c] = rst['comp1']

        if stop == 1:
            amp[:, c, :] = rst['amp']
        
        if cfg.prg == 1:
            print(f'Took {rst["iter"] + 1} iterations')

    comp = c_l.transpose((0, 2, 1))
    comp1 = c_sl.transpose((0, 2, 1))

    results.erp_new = 0
    results.residue = erp

    bl_wd = np.arange(-cfg.epoch_twd[0]/cfg.samp_interval,
                    -cfg.epoch_twd[0]/cfg.samp_interval+cfg.bl/cfg.samp_interval,
                    dtype=int)

    component = np.zeros((d1, d2, cfg.comp_num))
    component1 = np.zeros((d1, d2, cfg.comp_num))

    for j in np.arange(cfg.comp_num):
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
        results.residue = results.residue - component1[:, :, j]
        results.comps[cfg.comp_name[j]] = component[:, :, j]
        results.comps_sl[cfg.comp_name[j]] = component1[:, :, j]
        results.latencies[cfg.comp_name[j]] = cfg.comp_latency[j]*cfg.re_samp
        results.amps[cfg.comp_name[j]] = amp[:, :, j]

    results.comps[cfg.comp_name[0]] = baseline(results.comps[cfg.comp_name[0]],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)
    results.comps_sl[cfg.comp_name[0]] = baseline(results.comps_sl[cfg.comp_name[0]],bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)

    results.comps[cfg.comp_name[rst['trend_c']]] = baseline(results.comps[cfg.comp_name[rst['trend_c']]] + results.residue,bl_wd)
    results.comps_sl[cfg.comp_name[rst['trend_c']]] = baseline(results.comps_sl[cfg.comp_name[rst['trend_c']]] + results.residue,bl_wd)

    if  cfg.comp_num == 1:
        bl_wd = np.arange(-cfg.epoch_twd[0]/cfg.samp_interval, dtype='int')
        results.comps[cfg.comp_name[0]] = baseline(results.comps[cfg.comp_name[0]] + results.residue,bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)
        results.comps_sl[cfg.comp_name[0]] = baseline(results.comps_sl[cfg.comp_name[0]] + results.residue,bl_wd) + np.repeat(np.mean(erp[bl_wd, :], axis=0)[np.newaxis, :], epoch_length, axis=0)

    for j in np.arange(cfg.comp_num):
        results.erp_new = results.erp_new + results.comps[cfg.comp_name[j]]

    results.cfg = cfg0.copy()

    return results
