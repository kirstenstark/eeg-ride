function temp = move3(data,latency)
% for j = 1:size(data,3)
%     temp(:,:,j) = move2(data(:,:,j),round(latency(j)),'1d');
% end
% f = temp;
latency = latency(:)';
latency = round(latency);
[d1,d2,d3] = size(data);
temp = zeros(d1,d2,d3);


left = latency+1;left(left<=0) = 1;
right = d1+latency;right(right>d1) = d1;

left1 = -latency+1;left1(left1<=0) = 1;
right1 = d1-latency;right1(right1>d1) = d1;


for j = find(latency>-d1&latency<=d1)
    temp(left(j):right(j),:,j) = data(left1(j):right1(j),:,j);
end
