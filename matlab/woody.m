function [latency] = woody(data,cfg,dur)

[d1,d2,d3]=size(data);
latency = zeros(d3,1);
high_cutoff = fix(cfg.high_cutoff*10*d1*cfg.re_samp/1000);
bd = cfg.bd;
%low pass filter the data only for calculation of cross-correlation curve
for j = 1:d3    data(:,:,j) = filtering20(data(:,:,j),1,high_cutoff);end
%detrend data to remove drifting (only for latency estimation)
for j = 1:d3    data(:,:,j) = RIDE_detrend(data(:,:,j),[1,fix(d1*bd),fix(d1*(1-bd)),d1]);end
%remove the boundary effect
data = data.*repmat(RIDE_tukey(d1,bd*2),[1,d2,d3]);



data1 = move3(data,-latency);%figure;subplot(1,3,1);imagesc(mean(data1,3));
% try subplot(1,3,2);[a,b] = princomp(mean(data1,3));topoplot(a(:,1),chanlocs);end

% figure;imagesc(data(:,:,21));

if cfg.prg == 1 fprintf('estimating initial latency of C');end
clear temp5;clear temp11;
for iter = 1:1% Woody's method is only iterated once to 1. avoid arbitrary convergence to slow wave noise 2. comply with template approach
    if cfg.prg == 1 fprintf('.');end
    data1 = move3(data,-latency);
    
    if isfield(cfg,'template')
        if strcmp(cfg.template.method,'hanning')
        template = repmat(RIDE_hann(d1),[1,d2])...
            .*repmat(mean(mean(data,3)),[d1,1]);
        end
        if isfield(cfg.template,'hann_amp')
        template = repmat(RIDE_hann(d1),[1,d2])...
            .*repmat(cfg.template.hann_amp,[d1,1]);
        end
        if strcmp(cfg.template.method,'g_mean')
            template = cfg.temp;
        end
    end
    temp11 = [];
    
for m = 1:d3
    temp1 = (mean(data1(:,:,[1:d3]~=m),3));
    if isfield(cfg,'template') 
        if strcmp(cfg.template.method,'woody')
        template = (mean(data1(:,:,[1:d3]~=m),3));
        end
        temp1 = template;
    end
    
    temp = (data(:,:,m)); 
   
    for c = 1:d2
        temp11(:,c) = xcov(temp(:,c),temp1(:,c),fix(size(temp,1)/2),cfg.xc);%cross correlation
    end
    tem0 = mean(temp11(fix(size(temp11,1)/2-dur)+1:fix(size(temp11,1)/2+dur),:));

    temp5(:,m) = mean(temp11,2);

    temp5(:,m) = filtering10(temp5(:,m),1,fix(high_cutoff));
end

%make sure the highest point is not detected at the boundary
temp5 = RIDE_detrend(temp5,[1,fix(size(temp5,1)*bd),fix(size(temp5,1)*(1-bd)),size(temp5,1)]);
%!!!never apply time window here temp5 = temp5.*repmat(RIDE_tukey(size(temp5,1),bd*2),[1,size(temp5,2)]);


% disp(size(temp5));disp(dur);
temp5 = temp5(fix(size(temp5,1)/2-dur)+1:fix(size(temp5,1)/2+dur),:);


for j = 1:d3 
    if strcmpi(cfg.latency_search,'most_prob')
        latency(j) = nearest_latency(temp5(:,j), fix(size(temp5,1)/2));
    end
    if strcmpi(cfg.latency_search,'all')
        latency(j) = find_peak(temp5(:,j));
    end
end


% randomly assign the latency of the trials without x-corr peak
temp = ones(d3,1);
for j = 1:d3 
    if no_peak(temp5(:,j))==0 temp(j) = 0;end
end
latency(temp==0) = randn(length(find(temp==0)),1)*std(latency(temp==1)) + fix(size(temp5,1)/2);
latency(latency<1) = 1;latency(latency>d1)=d1;%make sure not exceed the boundary

latency = round(latency-median(latency));

% latency = find_peak2(temp5);
% latency = round(latency-median(latency));


% figure;imagesc(normalize(temp5));hold on;plot(latency + fix(size(temp5,1)/2),'o-');

% figure;imagesc(normalize(temp5));
end

% latency = find_peak2(squeeze(mean(data,2)));latency = round(latency-mean(latency));


if cfg.prg == 1 fprintf('done\n');end