function [Bmax, pMax] = Calc_Bmax(data, seg, K, theta, t, th, fastEst)
%CALC_BMAX calculates the highest segmentation benefit (B) w.r.t. the target's power (p).
%
%Description: 
% 	 This function finds the optimal segmentation benefit (B) in terms of target's power (p).
%    Details are introduced in CALC_PMAX.
% 
%Inputs: 
% 	 data (or y) - The data-cube to process (of the residual noise).
%    seg - Segmentation map.
%    K - Total scaling of the two segments. Null is 1.
%    theta - Total rotation angle of the two segments, in degrees. Null is 0.
% 	 t - A vector representing the spectrum of the target of interest.
%    th - The decision thresholds in terms of Pfa.
%    fastEst - A fast rough approximation by p1*sqrt(2).
% 
%Outputs: 
% 	 BMax - The segmentation benefit at the optimal target's power.
% 	 pMax - The optimal target's power.
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    [Bmax, pMax] = Calc_Bmax(data, seg, 4, 20, t, 0.01, 0);  % test spiting into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    fprintf('The optimal target's power is p_max=%g, and the respective segmentation benefit is B_max=%g\n', pMax, Bmax)
% 
%See also CALC_PMAX, CALC_B
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
if (nargin<7)
	fastEst = 0;
end

% calculates maximum of B(p)
[pMax, Bmax] = Calc_pmax(data, seg, K, theta, t, th, fastEst);
