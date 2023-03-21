function f = move2(aa,b,c)



b = round(b);
if nargin == 3
if strcmp(c,'1d') == 1 b = b*ones(1,size(aa,2));end
end


for j = 1:size(aa,2)
    aa(:,j) = move(aa(:,j),b(j),'move');
end


f = aa;


