import numpy as np


class RideCfg:

    def __init__(self,
                 comp_name,
                 comp_twd,
                 comp_latency,
                 sfreq,
                 epoch_twd=None,
                 rwd=200,
                 re_samp=2,
                 bd=0.2,
                 prg=1,
                 bl=200):

        assert len(comp_name) == len(comp_twd) == len(comp_latency), \
            '`comp_name`, `comp_twd`, and `comp_latency` must have the same length'

        comp_num = len(comp_name)

        samp_interval = 1000 / sfreq

        if epoch_twd is None:
            epoch_twd = np.array([-100, 1198])

        self.comp_name = comp_name
        self.comp_twd = comp_twd
        self.comp_latency = comp_latency
        self.comp_num = comp_num
        self.sfreq = sfreq
        self.samp_interval = samp_interval
        self.epoch_twd = epoch_twd
        self.rwd = rwd
        self.re_samp = re_samp
        self.bd = bd
        self.prg = prg
        self.bl = bl


    def copy(self):

        return RideCfg(comp_name=self.comp_name,
                       comp_twd=self.comp_twd,
                       comp_latency=self.comp_latency,
                       sfreq=self.sfreq,
                       epoch_twd=self.epoch_twd,
                       rwd=self.rwd,
                       re_samp=self.re_samp,
                       bd=self.bd,
                       prg=self.prg,
                       bl=self.bl)
