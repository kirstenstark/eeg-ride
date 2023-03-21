function f = interp2d(data,x1,x2,c)
for j = 1:size(data,2)
    temp(:,j) = interp1(x1,data(:,j),x2,c);
end
f = temp;