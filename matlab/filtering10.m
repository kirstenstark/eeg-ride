function f = filtering10(x,a,b)

for copyright = 1:1
% Code Author: Guang Ouyang, HKBU, 2010,
% The whole RIDE method developer: Guang Ouyang, Werner Sommer, Changsong Zhou (alphabetical order)
% Copyright (C) <2010>  Guang Ouyang, Werner Sommer, Changsong Zhou
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
end

%filter the data using ifft
%filter the data at a higher resolution (10 times as 'filtering')
%e.g., data = filtering10(data,23,35); 
%bandpass the data at [2.3Hz,3.5Hz] (assuming the data length is one second)









x = x(:);
x0 = mean(x);
x = x - x0;
n = 10*length(x);
Y = fft(x,n);
if b == 0
    H = [Y(1)',zeros(1,n-1)]';
end
if a == 0 && b~=0
    H = [Y(1:b+1)',zeros(1,n-(b)*2-3), Y(n-b-1:n)']';
end
if a ~= 0
    H = [zeros(1,a),Y(a+1:b+1)',zeros(1,n-b*2-1), Y(n-b+1:n-a+1)',zeros(1,a-1)]';
end
    f = real(ifft(H));
    f = f(1:length(x));
    
    if a==0 f = f+x0;end
