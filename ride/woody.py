import numpy as np
from .filtering import filtering20, filtering10
from .detrend import ride_detrend
from .tukey import ride_tukey
from .correct import move3
from .hann import ride_hann
from .nearest_latency import nearest_latency
from .find_peak import find_peak
from .no_peak import no_peak

def woody(data, cfg, dur):
    d1, d2, d3 = data.shape
    latency = np.zeros(d3)
    high_cutoff = int(cfg.high_cutoff * 10 * d1 * cfg.re_samp / 1000)
    bd = cfg.bd
    
    # Low-pass filter the data for cross-correlation calculation
    for j in range(d3):
        data[:, :, j] = filtering20(data[:, :, j], 1, high_cutoff)
    
    # Detrend data to remove drifting (only for latency estimation)
    for j in range(d3):
        data[:, :, j] = ride_detrend(
            data = data[:, :, j],
            twd = [0,
                   int(d1 * bd) - 1,
                   int(d1 * (1 - bd)) - 1,
                   d1 - 1])
    
    # Remove the boundary effect
    data *= np.tile(ride_tukey(d1, bd * 2).reshape(-1, 1, 1), (1, d2, d3))
    
    data1 = move3(data, -latency)
    
    if getattr(cfg, 'prg', 0) == 1:
        print('Estimating initial latency of C', end='')
    
    for iter in range(1):  # Woody's method is only iterated once
        if getattr(cfg, 'prg', 0) == 1:
            print('.', end='')
        
        data1 = move3(data, -latency)
        
        if hasattr(cfg, 'template'):
            if cfg.template.method == 'hanning':
                template = np.tile(ride_hann(d1).reshape(-1, 1), (1, d2)) * np.mean(np.mean(data, axis=2), axis=0)
            if hasattr(cfg.template, 'hann_amp'):
                template = np.tile(ride_hann(d1).reshape(-1, 1), (1, d2)) * cfg.template.hann_amp
            if cfg.template.method == 'g_mean':
                template = cfg.temp
        
        temp5 = np.zeros((d1, d3))
        for m in range(d3):
            temp1 = np.mean(data1[:, :, [i for i in range(d3) if i != m]], axis=2)
            if hasattr(cfg, 'template') and cfg.template.method == 'woody':
                template = temp1
                temp1 = template
            temp = data[:, :, m]
            
            temp11 = np.zeros((2 * (d1 // 2) + 1, d2))
            for c in range(d2):
                temp11[:, c] = xcov(temp[:, c],
                                    temp1[:, c],
                                    maxlags=int(np.fix(temp[:, c].shape[0]/2)),
                                    scaleopt=cfg.xc)
            
            temp5[:, m] = np.mean(temp11, axis=1)
            temp5[:, m] = filtering10(temp5[:, m], 1, high_cutoff)
        
        # Ensure highest point is not detected at the boundary
        temp5 = ride_detrend(temp5, [0,
                                     int(temp5.shape[0] * bd) - 1,
                                     int(temp5.shape[0] * (1 - bd)) - 1,
                                     temp5.shape[0] - 1])
        temp5 = temp5[d1 // 2 - dur : d1 // 2 + dur, :]
        
        for j in range(d3):
            if cfg.latency_search.lower() == 'most_prob':
                latency[j] = nearest_latency(temp5[:, j], temp5.shape[0] // 2)
            if cfg.latency_search.lower() == 'all':
                latency[j] = find_peak(temp5[:, j])
        
        temp = np.ones(d3)
        for j in range(d3):
            if no_peak(temp5[:, j]) == 0:
                temp[j] = 0
        
        latency[temp == 0] = np.random.randn(np.sum(temp == 0)) * np.std(latency[temp == 1]) + temp5.shape[0] // 2
        latency[latency < 1] = 1
        latency[latency > d1] = d1
        latency = round_like_matlab(latency - np.median(latency))
    
    if getattr(cfg, 'prg', 0) == 1:
        print('done')
    
    return latency
