function [a,b] = lls(x,y)


a = (mean(y)*sum(x.^2) - mean(x)*sum(x.*y))/(sum(x.^2)-length(x)*((mean(x)).^2));

b = (sum(x.*y) - length(x)*mean(x)*mean(y))/(sum(x.^2)-length(x)*((mean(x)).^2));