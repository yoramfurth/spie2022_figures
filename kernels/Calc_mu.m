function [mu, muAll, th2etaFun] = Calc_mu(data, seg, K, theta, t, p, MFtype)
%CALC_MU calculates the expectation of the with-target NSMF from zero, that is, given simulation generating factors only.
%
%Description: 
%    This function calculates the expectation from zero by computing all modules up to the required point. 
%    Note: if K or theta are not null, data and seg should contain two equal stationary segments, such as SIC data.
% 
%Inputs: 
% 	 data (or y) - The data-cube to process (of the residual noise).
%    seg - Segmentation map.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 	 t - A vector representing the spectrum of the target of interest. 
%    p - The target's power. 
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 mu - A single value representing the NSMF expectation.
% 	 muAll - NSMF expectation per segment ("Local" mode only).
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    [mu0, mu12] = Calc_mu(data, seg, 4, 20, t, 1, 'Local');  % test spliting into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    fprintf('The total expectation is - %g\n', mu0)
%    fprintf('The local expectations are - mu1=%g, mu2=%g\n', mu12(1), mu12(2))
% 
%See also EVAL_MU
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% calculates NSMF scores
y = Calc_y(data, seg, 1, 1);  % extracts the residual noise (if needed)
y = Split_y(y, seg, K, theta);  % splits into 2 distinct segments (if requested)
[MF_N, MF_T] = Calc_MF(y, t, MFtype, seg, p);

% calculates mu
[mu, muAll] = Eval_mu(MF_T, seg, MFtype);

% optionally creates handle to a th-to-eta inverse function
if (nargout>=3)
    th2etaFun = @(th) th2eta(MF_N(:),th);
end
