function y=distortion(x, gain, mix)
% Author: Dutilleux, Zölzer
%--------------------------------------------------------------------------
% This source code is provided without any warranties as published in 
% DAFX book 2nd edition, copyright Wiley & Sons 2011, available at 
% http://www.dafx.de. It may be used for educational purposes and not 
% for commercial applications without further permission.
%--------------------------------------------------------------------------
% Adapted by DTM 25/11/2020
% "Overdrive" simulation with symmetrical clipping
% x    - input
% gain - amount of distortion, >0
% mix  - mix of original and distorted sound, 1=only distorted
% Exam no: Y3858230

q=x*gain;
z=sign(q).*(1-exp(-abs(q)));
y=mix*z+(1-mix)*x;