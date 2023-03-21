function f = ms_rate(data,x,bin,y,c,varargin)
noplot = 0;
for j = 1:2:length(varargin)
    parameter = varargin{j};value = varargin{j+1};
    parameter = lower(parameter);
    switch parameter
        case 'bandwidth' 
            if ~isempty(value) bw = value;end
        case 'noplot'
            if ~isempty(value) noplot = value;end
    end
end %input parameters

if exist('bw','var')
    temp = ksdensity(data(~isnan(data)),x,'width',bw);
    f = temp*1000;
end


if ~exist('bw','var')
data = data(:);
temp = zeros(1,length(x));
for j = 1:length(data)
    for jj = 1:length(x)
        if ~isnan(data(j)) && data(j)>x(jj)-bin && data(j)<x(jj)+bin
            temp(jj) = temp(jj) + 1;
        end
    end
end
f = temp*1158/(2*bin*y);
end

if noplot == 0 plot(x,f,c);end
