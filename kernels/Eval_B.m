function B = Eval_B(AG, AL, AGnom, ALnom, dAGnom, dALnom, epsAG, epsAL)
%EVAL_B evaluates the segmentation benefit (B) 
%
%Description: 
%    Where the benefit is defined simply by B=AL/AG, some extreme cases require numeric approach. 
%    The first problem is on very low "th" where the denominator of "A" is too low. This is dealt 
%    by substituting AL,AG by their numerators. The second problem is very low AL,AG, which is 
%    dealt using L'Hôpital's approximation. See the thesis report Section 6.3. 
% 
%Inputs: 
% 	 AG, AL - Detection algorithm evaluator in global and local domain, respectively.
% 	 AGnom, ALnom - The numerator of AG and AL, respectively.
% 	 dAGnom, dALnom - Derivative of AGnom and ALnom, respectively. Usable for L'hopital's approximation.
%    epsAG, epsAL - Optional. What number is too low for AG and AL, respectively.
% 
%Outputs: 
% 	 B - The segmentation benefit. 
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
if (~exist('epsAG','var'))
    epsAG = 1e-5;
end
if (~exist('epsAL','var'))
    epsAL = 1e-3;
end

% numeric solution for approximating B=AL/AG
B = ALnom./AGnom;

% on very low values use L'Hôpital's approximation
k = (AG<epsAG & AL<epsAL);
if (any(k))
    B(k) = dALnom(k)./dAGnom(k);
end
