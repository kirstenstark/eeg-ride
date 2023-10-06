import numpy as np

from .detrend import ride_detrend
from .filtering import filtering20
from .helpers import mean_nan
from .median2d import median_2d
from .tukey import ride_tukey


def ride_iter(data, cfg):

    # # For debugging
    # cfg = cfg1
    # data = data[:, 61, :]

    d1, d2 = data.shape
    bd = cfg.bd

    data0 = data.copy()

    data = filtering20(data, 1, np.fix(10 * 20 * cfg.re_samp * d1 / 1000))
    # TODO: Compare with MNE filter

    max_latency = np.zeros((cfg.comp_num), dtype=int)
    min_latency = np.zeros((cfg.comp_num), dtype=int)
    length_c = np.zeros((cfg.comp_num), dtype=int)
    for j in range(cfg.comp_num):
        max_latency[j] = cfg.comp_latency[j].max()
        min_latency[j] = cfg.comp_latency[j].min()
        length_c[j] = d1 + max_latency[j] - min_latency[j]

    com_c = np.zeros((d1, cfg.comp_num))
    com_c1 = np.zeros((d1, d2, cfg.comp_num))
    amp_c = np.zeros((d2, cfg.comp_num))

    trend_i = cfg.comp_num-1
    stream_flow = np.arange(cfg.comp_num)
    for c in range(cfg.comp_num):
        if cfg.comp_name[c] == 'r':
            trend_i = c-1
    if trend_i == cfg.comp_num:
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i)])
    if trend_i == cfg.comp_num-2 : 
        stream_flow = np.concatenate([np.array([trend_i]), np.arange(trend_i), np.array([cfg.comp_num-1])])
    if trend_i == 0 : 
        stream_flow = np.array([1,0])

    l1 = np.zeros((cfg.inner_iter, cfg.comp_num))
    stop = 0
    stop_c = np.zeros((cfg.comp_num, 1))
    com_old = com_c.copy()

    # for iter in np.arange(1, 2): # For testing purposes, only runs iter = 1
    for iter in np.arange(cfg.inner_iter):

        if iter + 1 == cfg.inner_iter:
            stop = 1

        # track the convergence

        # !!!! CONTINUE HERE, AFTER FINISHING FIRST ITERATION OF FOR LOOP
        # TO DO: add if-loop for iter > 1

        if iter > 1: 
            for c in np.arange(cfg.comp_num):
                l1[iter-2,c]=np.sum(np.abs(com_c[:,c] - com_old[:,c]))
                if iter > 2:
                    if l1[iter-3,c] - l1[iter-2,c] < 0.01*(l1[0,c] - l1[1,c]): # convergence is now defined as 0.01
                        stop_c[c] = 1
                    if l1[iter-3,c] == l1[iter-2,c]:
                        stop_c[c] = 1
            if cfg.comp_num == 1:
                stop_c[c] = 1
        for c in np.arange(cfg.comp_num):
            com_old[:,c] = com_c[:,c]
            # decision of the termination of the interation of each component
            

        
        #for c in np.arange(cfg.comp_num):
            #com_old[:,c]=com_c[:,c] # decision of the termination of iteration of each component
        ## CAUTION: for-loop seems unnecessary

        # com_old = com_c.copy()

        for c in stream_flow:
            if stop_c[c] == 0:
                temp = data.copy()
                for j in np.arange(cfg.comp_num):
                    if j != c:
                        temp = temp-com_c1[:, :, j]
                residue = temp.copy()
                temp = np.empty((length_c[c], d2))
                temp[:] = np.nan
                for j in np.arange(d2):
                    temp[np.arange(-cfg.comp_latency[c][j]+max_latency[c],d1-cfg.comp_latency[c][j]+max_latency[c]),j] = residue[:,j]

                temp = ride_detrend(
                    temp,
                    np.array([0,
                              np.fix((cfg.comp_twd[c][1] - cfg.comp_twd[c][0]) * bd),
                              np.fix((cfg.comp_twd[c][1] - cfg.comp_twd[c][0]) * (1 - bd)) - 1,
                              cfg.comp_twd[c][1] - cfg.comp_twd[c][0] - 1], dtype=int) \
                                + max_latency[c] + int(cfg.comp_twd[c][0]))
                
                temp0 = median_2d(temp.T)
                temp0[np.concatenate([np.arange(cfg.comp_twd[c][0], dtype=int) + max_latency[c],
                                      np.arange(cfg.comp_twd[c][1] + max_latency[c]-1, length_c[c], dtype=int)])]=0
                temp0[np.arange(cfg.comp_twd[c][0]+ max_latency[c]-1, cfg.comp_twd[c][1] + max_latency[c], dtype=int)] = \
                    temp0[np.arange(cfg.comp_twd[c][0]+ max_latency[c]-1, cfg.comp_twd[c][1] + max_latency[c], 
                                    dtype=int)] * ride_tukey(cfg.comp_twd[c][1] - cfg.comp_twd[c][0] + 1, bd*2)
     
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
                for j in np.arange(cfg.comp_num):
                    if j != c:
                        temp = temp - com_c1[:,:,j]

                residue = temp.copy()
                temp = np.empty((length_c[c], d2))
                temp[:] = np.nan

                for j in np.arange(d2):
                    temp[np.arange(-cfg.comp_latency[c][j]+max_latency[c],
                                   d1-cfg.comp_latency[c][j]+max_latency[c]),
                         j] = residue[:,j]

                if cfg.final == 1:
                    temp[np.isnan(temp)] = 0
                    amp_c[:,c] = np.mean(
                        np.repeat(com_c[np.arange(cfg.comp_twd[c][0] - 1, cfg.comp_twd[c][1], dtype=int), c, np.newaxis], d2, axis=1) * temp[max_latency[c] + np.arange(cfg.comp_twd[c][0] - 1, cfg.comp_twd[c][1], dtype=int), :], axis=0)

            break
        # end of inner iter loop
    
    # for last iteration
    if cfg.final == 1:
        # release time window function and detrending
        # allocate trend to the last C component
        for c in stream_flow:
            temp = data0.copy()
            for j in np.arange(cfg.comp_num):
                if j != c:
                    temp = temp - com_c1[:,:,j]
                    # TODO: Probably only for microsaccades, temp-com_ms1 if 'ms' and 'var' exist
            residue = temp.copy()
            temp = np.empty((length_c[c], d2))
            temp[:] = np.nan

            for j in np.arange(d2):
                temp[np.arange(-cfg.comp_latency[c][j]+max_latency[c],
                               d1-cfg.comp_latency[c][j]+max_latency[c]),
                               j] = residue[:,j]
            temp0 = mean_nan(temp,1)

            if c==stream_flow[0]:
                temp0[np.arange(cfg.comp_twd[c][0] + max_latency[c],dtype=int)]=0.0
                tem = ride_tukey(cfg.comp_twd[c][1] - cfg.comp_twd[c][0] + 1, bd*2)
                tem = tem[np.arange(np.fix(len(tem)/2).astype(int))]
                temp0[np.arange(cfg.comp_twd[c][0] + max_latency[c]-1, 
                                cfg.comp_twd[c][0] + max_latency[c] + len(tem)-1, dtype=int)] = \
                                    temp0[np.arange(cfg.comp_twd[c][0] + max_latency[c]-1,\
                                                    cfg.comp_twd[c][0] + max_latency[c] + len(tem)-1, dtype=int)] * tem
                
            temp1 = np.repeat(temp0[:,np.newaxis],d2, axis=1)
            temp1[np.isnan(temp)] = np.nan
            temp1 = np.reshape(temp1.T[~np.isnan(temp1.T)], (d1, d2), order='F')
            com_c[:,c] = temp0[np.arange(max_latency[c], max_latency[c]+d1)]
            com_c1[:,:,c] = temp1.copy()
        
        c = stream_flow[0]
        temp = data0.copy()
        for j in np.arange(cfg.comp_num):
            if j != c:
                temp = temp - com_c1[:,:,j]
                # TODO: Probably only for microsaccades, temp-com_ms1 if 'ms' and 'var' exist
        residue=temp.copy()
        temp = np.empty((length_c[c], d2))
        temp[:] = np.nan
        for j in np.arange(d2):
            temp[np.arange(-cfg.comp_latency[c][j]+max_latency[c],
                           d1-cfg.comp_latency[c][j]+max_latency[c]),
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
