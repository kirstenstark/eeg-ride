function f = grouping(signal,num)
signal1 = signal;
for j = 1:size(signal,2)
%     disp(j);
%     if j>size(signal,2) - num
%         signal1(:,j) = mean(signal(:,[j-(j-(size(signal,2)-num)):j+num-(j-(size(signal,2)-num))]),2);
%     else
%     signal1(:,j) = mean(signal(:,[j:j+num]),2);
%     end
    a = j - round(num*j/size(signal,2));b = j + num - round(num*j/size(signal,2));
    signal1(:,j) = mean(signal(:,[a:b]),2);
end
f = signal1;