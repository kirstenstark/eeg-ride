function f = median_nan(x)
x = x(~isnan(x));
f = median(x);