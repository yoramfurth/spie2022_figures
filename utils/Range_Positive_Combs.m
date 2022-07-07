function aRng = Range_Positive_Combs(x, y)
%RANGE_POSITIVE_COMBS Calculates the range of angles of the positive cone. 
%
%Description: 
%    The range of angles of the positive cone is the range where any combination 
%    of (x,y) gives a vector with positive components. This is done by 
%    solving analytically the boundary condition. 
%  
%    Assuming we have a vector that already satisfies x>0 and we make linear combinations 
%    with another vector y which is situated outside the first cone, and which is not 
%    necessarily orthogonal to x. Then, all possible directions along those combinations 
%    might be described using a rotation angle "a" by: t(a)=cos(a)*x+sin(a)*y.
%    This can be rewritten by: t(a)=r*sin(a+phi), while r^2=x^2+y^2, and phi=atan(x/y).
%    Constraining a positive value means that "a" should be within [-phi,180-phi]+360k.
%
%Inputs: 
% 	 x - A multi-dimensional vector.
% 	 y - Another multi-dimensional vector (not parallel). 
% 
%Outputs: 
% 	 aRng - Range of angles where t(a) is positive. 
% 
%Example: 
%    a_rng = Range_Positive_Combs(x, y);
%    ta_rng = cosd(a_rng) .* [1;0] + sind(a_rng) .* [0;1];  % the two vectors on the boundary points
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

rngFun = @(rng1,rng2) rng1 : abs(rng2-rng1) : rng2;
rngFunVec = @(rngVec) rngFun(max(rngVec(:,1)), min(rngVec(:,2)));

phi = atan2d(x,y);
aRng = rngFunVec([0, 180] - phi); 
if (isempty(aRng))
    aRng = rngFunVec([0, 180] + mod(-phi,360));
end
