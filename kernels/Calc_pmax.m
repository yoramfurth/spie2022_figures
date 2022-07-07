function [pMax, Bmax] = Calc_pmax(data, seg, K, theta, t, th, fastEst)
%CALC_PMAX calculates the target's power (p) that provides the highest segmentation benefit (B).
%
%Description: 
% 	 This function finds for what target's power (p) the segmentation benefit (B) gets optimal.
%    The full solution starts with computing p-anchors, which is any "p" that satisfies mu==eta,
%    either globally or locally. Their bounds determines the range where B might be maximal. 
%    This range is then scanned and the maximum is estimated. CALC_PMAX has also a much faster 
%    approximation. Based on experimental observations the maximal power is estimated by just 
%    p1*sqrt(2). Although this is not accurate, its significant correlation to the real maximum 
%    makes it sufficient for a number of applications.  
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
% 	 pMax - The optimal target's power.
% 	 BMax - The segmentation benefit at the optimal target's power.
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    [pMax, Bmax] = Calc_pmax(data, seg, 4, 20, t, 0.01, 0);  % test spiting into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    fprintf('The optimal target's power is p_max=%g, and the respective segmentation benefit is B_max=%g\n', pMax, Bmax)
% 
%See also CALC_B
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

%% init.
if (nargin<7)
	fastEst = 0;
end

%% Finds the local p-anchors 
[~, qL, th2etaL] = Calc_mu(data, seg, K, theta, t, 1, 'Local');
pL = th2etaL(th) ./ qL;

if (fastEst)
	p1 = min(pL);  % p1 is the smallest local p-anchor
	pMax = p1 * sqrt(2);  % estimates p_max based on p1 anchor only
	if (nargout >= 2)
		Bmax = Calc_B(data, seg, K, theta, t, pMax, th);
	end
	return;
end

%% Finds the global p-anchor
[qG, ~, th2etaG] = Calc_mu(data, seg, K, theta, t, 1, 'Global');
pG = th2etaG(th) / qG;

%% Scans B(p) in the range of p-anchors (plus factor 10 margins)
pAnc = [pL(:); pG];  % all p-anchors
minRng = min(pAnc)/10;
maxRng = max(pAnc)*10;
pRng = 10.^(log10(minRng): 0.01 :log10(maxRng));
B = Calc_B(data, seg, K, theta, t, pRng, th);

%% Estimates the maximum of B(p)
B0 = imfilter(B, fspecial('gaussian',19,3), 'symmetric');
[B0max, kMax] = max(B0);
pMax = pRng(kMax);
Bmax = max(B0max, B(kMax));

