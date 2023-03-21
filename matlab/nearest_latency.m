function f = nearest_latency(signal,a)
% signal = signal + rand(size(signal,1),size(signal,2))/100000000;

% be careful about the 'peak' with equal adjacent values

temp = zeros(1,length(signal));
for j = 2:length(signal)-1
    if signal(j+1)<signal(j) && signal(j)>signal(j-1)
        temp(j) = 1;
    end
end
if sum(temp)==0 f = a;end
if sum(temp)~=0 
temp1 = find(temp==1);temp2 = abs(temp1-a);
f = temp1(temp2==min(temp2));
f = f(1);
end


% 
% if a == 1 && signal(2)>signal(1)
%     for j = 1:length(signal)-1
%         if signal(j+1)<signal(j) latency = j;break;end
%         if j == length(signal)-1 latency = j+1;end
%     end
% end
% if a == 1 && signal(2)<signal(1)
%     for j = 2:length(signal)-1
%         if signal(j+1)<signal(j)&&signal(j)>signal(j-1) latency = j;break;end
%         if j == length(signal)-1
%             if signal(1) > signal(j+1) latency = 1;end
%             if signal(1) < signal(j+1) latency = j;end
%         end
%     end
% end
% if a == length(signal)
%     a = 1;
%     signal = signal(abs([1:length(signal)]-length(signal)-1));
%     if a == 1 && signal(2)>signal(1)
%       for j = 1:length(signal)-1
%         if signal(j+1)<signal(j) latency = j;break;end
%         if j == length(signal)-1 latency = j+1;end
%       end
%     end
%         if a == 1 && signal(2)<signal(1)
%             for j = 2:length(signal)-1
%                 if signal(j+1)<signal(j)&&signal(j)>signal(j-1) latency = j;break;end
%                 if j == length(signal)-1
%                     if signal(1) > signal(j+1) latency = 1;end
%                     if signal(1) < signal(j+1) latency = j;end
%                 end
%             end
%         end
%         latency = length(signal) - latency + 1;
%         
% end
% 
% if a > 1 && a < length(signal)
%     if signal(a+1)<signal(a) && signal(a-1)<signal(a) latency = a;end
%     if signal(a+1)>signal(a)
%         jj = a;
%         while 1
%             jj = jj+1;
%             if jj == length(signal) latency = jj;break;end
%             if signal(jj) < signal(jj-1) latency = jj-1;break;end
%         end
%     end
%     if signal(a-1)>signal(a)
%         jj = a;
%         while 1
%             jj = jj-1;
%             if jj == 1 latency = jj;break;end
%             if signal(jj) < signal(jj+1) latency = jj + 1;break;end
%         end
%     end
% end
% f = latency;