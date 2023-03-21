function f = single_trial_RIDE(data,results,name,channel,varargin)

if strcmpi(name,'erp')
    f = squeeze(data(:,channel,:));
    return;
end


comp_name = results.cfg.comp.name;
temp = squeeze(data(:,channel,:));
si = results.cfg.samp_interval;

c_ind = find(strcmp(results.cfg.comp.name,name));
twd = fix((results.cfg.comp.twd{c_ind}-results.cfg.epoch_twd(1))/si);
if strcmp(name,'r')
    twd = results.cfg.comp.twd{end}+median(results.latency0{end});
    twd(twd<200) = 200;
    twd(twd>results.cfg.epoch_twd(2)) = results.cfg.epoch_twd(2);
    twd = fix((twd-results.cfg.epoch_twd(1))/si);
end
    


bd = results.cfg.bd;

for j = 1:size(temp,2)
    for k = 1:length(results.cfg.comp.name)
        if ~strcmp(results.cfg.comp.name{k},name)
            eval(['temp(:,j) = temp(:,j) - move(results.',comp_name{k},'(:,channel),fix(results.latency_',comp_name{k},'(j)/si),''move'');']);
        end
    end
end

temp1 = temp;
% disp(d);
for j = 1:length(varargin)
    if strcmp(varargin{1},'synced')
    eval(['temp1 = move2(temp,-fix((results.latency_',name,')/si),''2d'');']);

    end
    
    if strcmp(varargin{1},'refined')
    eval(['temp1 = move2(temp,-fix((results.latency_',name,')/si),''2d'');']);
    d = [twd(1) twd(1)+fix(twd(2)-twd(1))*bd twd(1)+fix(twd(2)-twd(1))*(1-bd) twd(2)];

    temp1 = RIDE_detrend(temp1,fix(d));

    temp1([1:twd(1),twd(2):size(temp1,1)],:) = 0;
    temp1(twd(1):twd(2),:) = temp1(twd(1):twd(2),:).*repmat(RIDE_tukey(twd(2)-twd(1)+1,bd*2),[1,size(temp1,2)]);
    eval(['temp1 = move2(temp1,fix((results.latency_',name,')/si),''2d'');']);
    end
    
    
    if strcmp(varargin{1},'refined_normalized')
    eval(['temp1 = move2(temp,-fix((results.latency_',name,')/si),''2d'');']);
    d = [twd(1) twd(1)+fix(twd(2)-twd(1))*bd twd(1)+fix(twd(2)-twd(1))*(1-bd) twd(2)];
    
    temp1 = normalize(temp1,twd(1):twd(2));
    temp1 = RIDE_detrend(temp1,fix(d));

    temp1([1:twd(1),twd(2):size(temp1,1)],:) = 0;
    temp1(twd(1):twd(2),:) = temp1(twd(1):twd(2),:).*repmat(RIDE_tukey(twd(2)-twd(1)+1,bd*2),[1,size(temp1,2)]);
    eval(['temp1 = move2(temp1,fix((results.latency_',name,')/si),''2d'');']);
    end
    
    if strcmp(varargin{1},'normalized')
    eval(['temp1 = move2(temp,-fix((results.latency_',name,')/si),''2d'');']);
    d = [twd(1) twd(1)+fix(twd(2)-twd(1))*bd twd(1)+fix(twd(2)-twd(1))*(1-bd) twd(2)];
    
    temp1 = normalize(temp1,twd(1):twd(2));
%     temp1 = RIDE_detrend(temp1,fix(d));

%     temp1([1:twd(1),twd(2):size(temp1,1)],:) = 0;
%     temp1(twd(1):twd(2),:) = temp1(twd(1):twd(2),:).*repmat(RIDE_tukey(twd(2)-twd(1)+1,bd*2),[1,size(temp1,2)]);
    eval(['temp1 = move2(temp1,fix((results.latency_',name,')/si),''2d'');']);
    end
    
    if strcmp(varargin{1},'syn_refined')
    eval(['temp1 = move2(temp,-fix((results.latency_',name,')/si),''2d'');']);

    d = [twd(1) twd(1)+fix(twd(2)-twd(1))*bd twd(1)+fix(twd(2)-twd(1))*(1-bd) twd(2)];

    temp1 = RIDE_detrend(temp1,fix(d));

    temp1([1:twd(1),twd(2):size(temp1,1)],:) = 0;
    temp1(twd(1):twd(2),:) = temp1(twd(1):twd(2),:).*repmat(RIDE_tukey(twd(2)-twd(1)+1,bd*2),[1,size(temp1,2)]);
    end
    
end

f = temp1;

