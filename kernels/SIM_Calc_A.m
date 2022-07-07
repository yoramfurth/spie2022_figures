function [A, Anom, AnomLop] = SIM_Calc_A(phiG, K, theta, t, th, MFtype)
%SIM_CALC_A calculates A(p) analytically given simulation generating factors.
%
%Description: 
% 	 The metric "A" evaluates how much relative superior detection is achieved over a range of thresholds.
%    This function calculates it analytically from zero, that is, given simulation generating factors only.
%    It does it by computing all modules up to the required point. 
% 
%Inputs: 
% 	 phiG - The global covariance matrix.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 	 t - A vector representing a target's spectrum. 
% 	 th - Decision threshold in terms of Pfa.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 A(p) - Evaluation of the detection algorithm as a function of p.
% 	 Anom(p) - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom(p) - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 
%Example: 
%    AR = 120; % major-to-minor aspect-ratio of the covariance matrix
%    phiG = SIM_Get_phi(AR);
%    t = [20 90]';
%    p = 0.02;  % target's power to test
%    th = 0.01;  % desired number of false alarms
%    A = SIM_Calc_A(phiG, 4, 20, t, th, 'Local');  % determines 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    fprintf('A(th)=%g for a threshold of th=%g, and a target's power of p=%g\n', A(p), th, p);
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

%step 1: splits data
[phi1, phi2] = SIM_Split_phi(phiG, K, theta);  % extract the data background

%step 2: calculates MF distributions
[fN,~,~,cN,cT,icN] = SIM_Calc_PDF(phi1, phi2, t, 0.5, MFtype);

%step 3: Calculates A(p)
eta = -icN(th);
[A, Anom, AnomLop] = SIM_Eval_A(cN, cT, fN, eta); 

