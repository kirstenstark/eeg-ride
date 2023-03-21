function f = no_peak(data)
temp=0;
for j = 2:length(data)-1
    if data(j-1)<data(j) && data(j)>data(j+1)
        temp = 1;
    end
end
f = temp;