function f = filtering20(x,a,b)
f = x;
for j = 1:size(x,2)
    f(:,j) = filtering10(x(:,j),a,b);
end
