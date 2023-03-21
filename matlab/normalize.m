function f = normalize(x,interval)


if isvector(x)
    f = (x - mean(x(~isnan(x))))/std(x(~isnan(x)));
end

if ~isvector(x)
temp = mean(x);
temp1 = std(x);
f = (x-temp(ones(1,size(x,1)),:,:))./temp1(ones(1,size(x,1)),:,:);
end




if nargin == 2
    temp = mean(x(interval,:,:));
temp1 = std(x(interval,:,:));
f = (x-temp(ones(1,size(x,1)),:,:))./temp1(ones(1,size(x,1)),:,:);
end