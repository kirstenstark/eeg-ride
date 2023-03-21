function f = reference(data)
d = size(data,2);
temp = mean(data,2);
f = data - temp(:,ones(1,d),:);        
        