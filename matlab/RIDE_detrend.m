function f = RIDE_detrend(data,twd)


% index = isnan(data);
% data(index) = 0;
% [d1,d2,d3] = size(data);
% for j = 1:d2
%     for k = 1:d3
%         a = mean(data(twd(1):twd(2),j,k));
%         b = mean(data(twd(3):twd(4),j,k));
%         temp = ([1:d1])*(b-a)/((twd(4)+twd(3))/2-(twd(2)+twd(1))/2);
%         temp = temp - mean(
%         data(:,j,k) = data(:,j,k) - temp';
%     end
% end
% f = data;
% f(index) = nan;





index = isnan(data);
data(index) = 0;
[d1,d2,d3] = size(data);
d = ((twd(4)+twd(3))/2-(twd(2)+twd(1))/2);

temp0 = repmat([1:d1]',[1,d2,d3]);

        a = mean(data(twd(1):twd(2),:,:),1);
        b = mean(data(twd(3):twd(4),:,:),1);
        temp = (b-a)/d;
        data = data - temp0.*temp(ones(1,d1),:,:);

        temp = mean(data(twd(1):twd(2),:,:),1);
        data = data - temp(ones(1,d1),:,:);
        
f = data;
f(index) = nan;