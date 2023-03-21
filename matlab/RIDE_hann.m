function f = RIDE_hann(d)
f = 0.5*(1-cos(2*pi*((1:d)-1)/(d-1)));
f = f(:);