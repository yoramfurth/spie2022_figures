function [tmax, KbMax] = SIM_Calc_tmax(phiG, phiS)
%SIM_CALC_TMAX calculates the optimal Kb(t) using whitening domain.
%
%Description: 
%    The function implement an estimates analytically the optimal target's direction, in terms of segmentation
%    benefit (B). The algorithm maximizes Kb(t), which incorporates the second moment of "B". In comes that  
%    all possible "Kb" values for one segment (e.g. "phiS{1}") are represented by the ellipsoid of "phiG" 
%    whitened  by this segment ("phiG_w"). The biggest "Kb" in this case is where the largest radius occurs, 
%    which is just on the major eigenvector of this whitened matrix ("EG_w(:,1)"). This process is repeated  
%    for every segment, and the highest major eigenvector is selected. Note that this technique ignore that 
%    target's power, implying that it refers to the best power on each direction. It practice this calculation    
%    is a super fast comparing to other alternatives, and it provides a pretty close estimation. More details  
%    are available in "readme.txt" and in our related paper. 
% 
%Inputs: 
% 	 phiG - Global covariance matrix.
% 	 phiS - Cell array of local covariance matrixes.
% 
%Outputs: 
% 	 tmax - The optimal target's direction, in terms of segmentation benefit.
% 	 KbMax - Directional inhomogeneity impact at "tmax".
% 
%Example: 
%    AR = 120; % major-to-minor aspect-ratio of the covariance matrix
%    phiG = SIM_Get_phi(AR);
%    [phi1, phi2] = SIM_Split_phi(phiG, 4, 20);  % determines 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    [tmax, KbMax] = SIM_Calc_tmax(phiG, {phi1, phi2})
%    fprintf('The highest Kb is %g, and it occurs at t=[%g,%g]\n', KbMax, tmax(1), tmax(2));
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Whitens the global covariance (phiG) matrix by each of the local matrixes (per phiS{})
sMax = numel(phiS);  % number of segments
d_w_max = zeros(1,sMax);
EG_w = cell(1,sMax); % Global eigenvectors, whitened
for ks = 1:sMax
    w = phiS{ks}^-0.5; % whitening operator
    phiG_w = w' * phiG * w;  % whitened phiG
    [EG_w{ks}, d_w] = eig_sorted(phiG_w);
    d_w_max(ks) = d_w(1);  % the maximal eigenvalue
end

% The global maximum 
[KbMax, sMax] = max(sqrt(d_w_max));  % this maximum among all the segments
tmax = phiS{sMax}^0.5 * EG_w{sMax}(:,1);  % inverse whitening. Reverts to the original domain.
