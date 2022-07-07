function [mu, muAll] = SIM_Calc_mu(phi1, phi2, t, p, NSeg, MFtype)
%SIM_CALC_MU Calculates analytically the expectation of the with-target NSMF.
%
%Description: 
%    This function provides one single value representing the NSMF expectation  
%    with-target. Based on EVAL_MU, it calculates the expectations directly
%    from NSMF components. The related background appears in EVAL_MU and is
%    detailed in readme.txt and in the related thesis report Section 5.4.2.
% 
%Inputs: 
% 	 phi1, phi2 - The covariance matrix of segment #1 and #2, respectively.
% 	 t - A matrix representing an array of vectors representing targets' spectrum. 
%    p - The power of the target. 
%    NSeg - Number of pixels per segment.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 mu - A single value representing the NSMF expectation.
% 	 muAll - NSMF expectation per segment ("Local" mode only).
% 
%Example: 
%    phi1 = [1 0; 0 0.25];
%    phi2 = [0.25 0; 0 1];
%    t = [20 90]';
%    [mu0, mu12] = SIM_Calc_mu(phi1, phi2, t, 0.02, 0.5, 'Local');
%    fprintf('The total expectation is - %g\n', mu0)
%    fprintf('The local expectations are - mu1=%g, mu2=%g\n', mu12(1), mu12(2))
% 
%See also EVAL_MU
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
NSeg = NSeg.*[1 1]; % number of pixels per segment
nSeg = NSeg/sum(NSeg); % apriori probability per segment

% selects between global and local mode
if (strcmpi(MFtype,'Global'))
    mu = SIM_Calc_mu_NGMF(phi1, phi2, t, p, nSeg);
    muAll = [];
	
else  % the local case - per segment calculation
    [mu, muAll] = SIM_Calc_mu_NSMF(phi1, phi2, t, p, nSeg);
end


function muG = SIM_Calc_mu_NGMF(phi1, phi2, t, p, nSeg)
%SIM_CALC_MU_NGMF calculates the expectation for the global case.

% calculates the global covariance matrix
phiG = nSeg(1)*phi1 + nSeg(2)*phi2;

% calculates the global expectation
muG = p * Eval_q(t, phiG);


function [mu0, muAll] = SIM_Calc_mu_NSMF(phi1, phi2, t, p, nSeg)
%SIM_CALC_MU_NSMF calculates a whole expectation for the local case.

% calculates expectation per segment
mu1 = p * Eval_q(t, phi1); 
mu2 = p * Eval_q(t, phi2); 
muAll = [mu1; mu2]; 

% calculates the neutral equivalence, i.e. the expectation in the neutral case where  
% all the expectations gets identical. (see the thesis report Section 5.4.2).
z_parallel = @(r) 1./sum(1./r);  % calculates parallel impedance
mu0 = sqrt(z_parallel(muAll(:).^2 ./ nSeg(:)));

