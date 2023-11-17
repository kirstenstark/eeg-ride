import matplotlib.pyplot as plt
import numpy as np


class RideResults:

    def __init__(self,
                 erp,
                 latency0,
                 erp_new=None,
                 residue=None,
                 comps=None,
                 comps_sl=None,
                 latencies=None,
                 amps=None,
                 cfg=None):

        if comps is None:
            comps = {}

        if comps_sl is None:
            comps_sl = {}

        if latencies is None:
            latencies = {}

        if amps is None:
            amps = {}

        self.erp = erp
        self.latency0 = latency0
        self.erp_new = erp_new
        self.residue = residue
        self.comps = comps
        self.comps_sl = comps_sl
        self.latencies = latencies
        self.amps = amps
        self.cfg = cfg

    def plot(self):
        """Plots time courses of the ERP and RIDE components."""

        n_comps = len(self.comps)
        ncols = 1 + n_comps
        fig, axs = plt.subplots(1, ncols, figsize=(10, 3))

        x = np.linspace(self.cfg.epoch_twd[0],
                        self.cfg.epoch_twd[1],
                        self.erp.shape[0])

        y_lim = self._get_ylim()

        axs[0].plot(x, self.erp, linewidth=0.7)
        axs[0].set_title('erp')
        axs[0].set_xlabel('Time')
        axs[0].set_ylabel('Amplitude')
        axs[0].set_ylim(y_lim)

        for ax, key in zip(axs[1:], self.comps):
            ax.plot(x, self.comps[key], linewidth=0.7)
            ax.set_title(key)
            ax.set_ylim(y_lim)
            ax.set_xlabel('Time')

        return fig

    def _get_ylim(self):
        """Gets appropriate y-axis limits for the time course plot."""

        y_max_erp = np.abs(self.erp.max())
        y_max_comps = [np.abs(comp.max()) for comp in self.comps.values()]
        y_max = max(y_max_erp, *y_max_comps)

        # Round up "small" numbers (e.g., volts) to nearest decimal
        if y_max < 1.0:
            n_digits = int(np.ceil(-np.log10(np.abs(y_max))))
            y_max = np.ceil(y_max * 10 ** n_digits) / 10 ** n_digits
        else:  # Round up "large" numbers (e.g., microvolts) to nearest integer
            y_max = np.ceil(y_max)

        return [-y_max, y_max]
