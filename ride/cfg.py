import numpy as np


def ride_cfg(comp,
             samp_interval,
             epoch_twd=np.array([-100, 1198]),
             rwd=200,
             re_samp=2,
             bd=0.2,
             prg=1,
             bl=200):

    # Continue making this a user-facing function

    cfg = {'samp_interval': samp_interval,
           'epoch_twd': epoch_twd,
           'comp': comp}
    cfg['comp_num'] = len(cfg['comp']['name'])
    cfg['rwd'] = rwd
    cfg['re_samp'] = re_samp
    cfg['bd'] = bd  # Alpha value for Tukey window
    cfg['prg'] = prg  # Print progress to console
    cfg['bl'] = bl  # Baseline length
