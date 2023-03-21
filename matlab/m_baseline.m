function f = m_baseline(x,interval)
index = isnan(x);
x(isnan(x))=0;

if nargin == 1
temp = median_3d(x);
temp = reshape(temp,1,size(temp,1),size(temp,2));
f = x-temp(ones(1,size(x,1)),:,:);
end

if nargin == 2
    
temp = median_3d(x(interval,:,:));
temp = reshape(temp,1,size(temp,1),size(temp,2));
f = x-temp(ones(1,size(x,1)),:,:);
end
f(index) = nan;