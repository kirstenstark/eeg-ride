
function y = median_2d(x)
  x(isnan(x)) = 0;%zeros padding
  
  
  
  
  
  
  
  
%   
%   
%   s = size(x);
%   % Sort along given dimension
%   x = sort(x);
%   total_dots = sum(ones(s));
% 
%       half = round(total_dots/2);
%       half = half+([1:s(2)]-1)*s(1);
%       y = x(half);
%       y = y(:);
% 
%   
%       
      
      
        
  
  s = size(x);
  % Sort along given dimension
  x = sort(x);
  
      y = x(round(s(1)/2) + ([1:s(2)]-1)*s(1));
      y = y(:);