function y=overdrive(x)
% Author: Dutilleux, Zölzer
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------
% Adapted by DTM 25/11/2020
% "Overdrive" simulation with symmetrical clipping
%  x - input
% Last edited on 31/12/2021
% Exam no: Y3858230

N=length(x);
th=1/3; % threshold for symmetrical soft clipping 

for i=1:N
   if abs(x(i))< th 
       y(i)=2*x(i);
   end
   if abs(x(i))>=th 
      if x(i)> 0 
          y(i)=(3-(2-x(i)*3).^2)/3;
      end
      if x(i)< 0
          y(i)=-(3-(2-abs(x(i))*3).^2)/3; 
      end
   end  
   if abs(x(i))>2*th 
      if x(i)> 0
          y(i)=1;
      end
      if x(i)< 0
          y(i)=-1;
      end
   end
end