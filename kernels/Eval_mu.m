function [mu, muAll] = Eval_mu(MF_T, seg, MFtype)
%EVAL_MU Evaluates the expectation of the with-target NSMF.
%
%Description: 
%    This function provides one single value representing the NSMF expectation  
%    with-target. Where in the global it is trivial, in the local case it requires
% 	 to decompose all the per-segment expectations. This function does it by:
%             mu0 = 1/sqrt(sum(nS/muS^2) ,
%    where "nS" is the relative size of segment number "s" and "muS" is the 
%    expectation within that segment, actually the average. The main justification 
%    is the inverse connection between any NSMF expectation and its related 
%    covariance matrix, given by mu=p*sqrt(t'*phi^-1*t), while mu is expectation  
%    and "phi" is its relate covariance matrix. Combining that phiG=sum(n_s*phiS), 
%    derives the above expression. see more in readme.txt and in the related thesis  
%    report Section 5.4.2.
% 
%Inputs: 
% 	 MF_T - The NSMF with target. 
%    seg - Segmentation map. Required in the "Local" mode only. 
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 mu - a single value representing the NSMF expectation.
% 	 muAll - NSMF expectation per segment ("Local" mode only).
% 
%Example: 
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    [NSMF_N, NSMF_T] = Calc_MF(y, t, 'Local', seg, 0.05);  % apply the NSMF
%    [mu0, mu12] = Eval_mu(NSMF_T, seg, 'Local');
%    fprintf('The total expectation is - %g\n', mu0)
%    fprintf('The local expectations are - mu1=%g, mu2=%g\n', mu12(1), mu12(2))
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
if (strcmpi(MFtype,'Global'))
	seg = [];  
end

% selects between global and local mode
if (strcmpi(MFtype,'Global') || isempty(seg))  % the global case
    mu = mean(MF_T);  % evaluate the global expectation by just averaging
    muAll = [];
    return 
    
else  % the local case - per segment calculation
    % estimates the MF expectation per segment
    sAll = unique(seg(:))';
    muAll = zeros(numel(sAll),1); 
    for s = sAll
        k = (seg==s);  % mask of the current segment
        MF_T_inSeg = MF_T(k(:));  % with-target MF of this segment's pixels 
        muAll(s,:) = mean(MF_T_inSeg);  % evaluate the expectation per-segment by averaging
    end
    
    % calculates the neutral equivalence, i.e. the expectation in the neutral case where  
    % all the expectations gets identical. (see the thesis report Section 5.4.2).
    NSeg = accumarray(seg(:),1);  % number of pixels per segment
    nSeg = NSeg / sum(NSeg);  % apriori probability per segment
	z_parallel = @(r) 1./sum(1./r);  % calculates parallel impedance
    mu0 = sqrt(z_parallel(muAll.^2 ./ nSeg));
	mu = mu0;
end
