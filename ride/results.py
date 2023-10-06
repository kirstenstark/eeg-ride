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
