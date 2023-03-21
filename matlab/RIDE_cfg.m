function cfg = RIDE_cfg(cfg)

if ~isfield(cfg,'ave_refer') 
    cfg.ave_refer = 0;%average reference, default: no
end
if ~isfield(cfg,'re_samp') 
    cfg.re_samp = cfg.samp_interval;%re-sampling, default: no
end
if ~isfield(cfg,'high_cutoff') 
    cfg.high_cutoff = 4;%high cutoff for cross-correlation curve, default: 5
end
if ~isfield(cfg,'bd') 
    cfg.bd = 0.2;%alpha value for tukey window,also bd/2 is the length of edge for detrending
end
if ~isfield(cfg,'bl') 
    cfg.bl = 200;%baseline time window
end
if ~isfield(cfg,'rwd') 
    cfg.rwd = 200;%minimal left boundary of R time window
end
cfg.comp_num = length(cfg.comp.name);
if ~isfield(cfg,'xc')
    cfg.xc = 'coeff';%cross-covariance
end
if isfield(cfg,'template')
    if ~isfield(cfg.template,'method')
        cfg.template.method = 'woody';
    end
end
if ~isfield(cfg,'latency_search') 
    cfg.latency_search = 'most_prob';%average reference, default: no
end

if ~isfield(cfg,'prg') 
    cfg.prg = 1;%show progress
end

    