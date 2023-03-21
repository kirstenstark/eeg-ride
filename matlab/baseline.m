function f = baseline(x,a)

x(isnan(x))=0;

if nargin == 1
temp = mean(x);
f = x-temp(ones(1,size(x,1)),:,:);
end


if nargin == 2
%     for j = 1:size(x,2)
%         for k = 1:size(x,3)
%             f(:,j,k) = x(:,j,k) - mean(x(~isnan(x(:,j,k)),j,k));
%         end
%     end
    temp = mean(x(a,:,:));
    f = x-temp(ones(1,size(x,1)),:,:);
end



