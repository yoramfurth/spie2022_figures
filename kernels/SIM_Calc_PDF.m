function [fN,fT,z,cN,cT,icN] = SIM_Calc_PDF(phi1, phi2, t, NSeg, MFtype)
%SIM_CALC_PDF calculates decision distributions of NSMF along the discriminant axis.
%
%Description: 
% 	 While applying NSMF on data without target the scores have one distribution (fN),  
%    and when a target is appended they have another distribution (fT). This function 
%    calculates distributions that correspond to the second moment only, represented by
%    the covariance matrixes of each segment (phi1, phi2). In the global domain this
%    gives a zero mean GMM for fN, and a equal-shape shifted GMM for fT. In the local 
%    domain this give a standard Gaussian for fN, and a dual mode GMM for fT. See more
%    details in the thesis report Section 2.4.
%
%    In order to allow fast analytic calculations, this function returns handles to  
%    functions, rather than concrete numbers. Their arguments are the scores along 
%    the discriminant axis (z), and the target's power (p). It additionally outputs also
%    an inverse tranform function, to transfer from the Pfa domain (given by cumulative
%    false alarms, and denoted "th") back to the scores domain. Note that the code 
%    below denotes "m" for mean and "s" for std.
% 
%Inputs: 
% 	 phi1, phi2 - The covariance matrix of segment #1 and #2, respectively.
% 	 t - A matrix representing an array of vectors representing targets' spectrum. 
%    NSeg - Number of pixels per segment.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 fN(z) - PDF of the scores of NSMF with no target.
% 	 fT(z,p) - PDF of the scores of NSMF with target.
% 	 z(p) - Ticks along the discriminant axis.
% 	 cN(z) - CDF of the scores of NSMF with no target.
% 	 cT(z,p) - CDF of the scores of NSMF with target.
% 	 icN(th) - inverse CDF function of the scores of NSMF with no target.
% 
%Example: 
%    phi1 = [1 0; 0 0.25];
%    phi2 = [0.25 0; 0 1];
%    t = [20 90]';
%    [fNz, fTz, z] = SIM_Calc_PDF(phi1, phi2, t, 0.5, 'Local'); 
%    p = 0.02;  % target's power to test
%    th = 0.01;  % desired number of false alarms
%    z = z(p); 
%    figure; plot(z,[fN(z); fTz(z,p)],'b');
%    legend('without target','with target');
% 
%See also NORMPDF
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
NSeg = NSeg.*[1 1]; % size of population of each segment
nSeg = NSeg/sum(NSeg); % apriori probabilities each segment

% selects between global and local mode
if (strcmpi(MFtype,'Global'))
    [fN,fT,z,cN,cT,icN] = SIM_Calc_PDF_NGMF(phi1, phi2, t, nSeg);

else  % the local case - per segment calculation
    [fN,fT,z,cN,cT,icN] = SIM_Calc_PDF_NSMF(phi1, phi2, t, nSeg);
end


function [fN,fT,z,cN,cT,icN] = SIM_Calc_PDF_NGMF(phi1, phi2, t, nSeg)
%SIM_CALC_PDF_NGMF calculates decision distributions for the global case.
% In the global case we have difference variances while evaluating per segment, but the the same shape for w\wo target.

% Determines GMM components
phiG = nSeg(1)*phi1 + nSeg(2)*phi2;
q = Eval_q(t, phiG); % the expectation for a full target  
MF_STD = q;  % by luck, in the NSMF algorithm, the STD of the matched filter is equal to the target's expectation, so we can reuse it
s(1) = sqrt(t'/phiG * phi1 * (phiG\t)) / MF_STD; 
s(2) = sqrt(t'/phiG * phi2 * (phiG\t)) / MF_STD;

% Generates samples 
z = GenerateSamples(q, s);  % z is a function, z(p) 

% Builds GMM
s_1 = 1./s;  % "_1" denotes "1/"
fN = @(z) (nSeg.*s_1) * normpdf(z .* s_1');
cN = @(z) nSeg * normcdf(z .* s_1');
fT = @(z,p) (nSeg.*s_1) * normpdf((z - p*q) .* s_1');
cT = @(z,p) nSeg * normcdf((z - p*q) .* s_1');

% Builds "th" inverse function, i.e. eta(th)
icN = @(th) fzero(@(eta)cN(eta) - th, 0);
icN = @(th) fevalInRange(icN, th, 0, 1-eps);


function [fN,fT,z,cN,cT,icN] = SIM_Calc_PDF_NSMF(phi1, phi2, t, nSeg)
%SIM_CALC_PDF_NSMF calculates decision distributions for the local case (i.e. per segment).
% In the local case we have variance 1 for any segment regardless w\wo target, whereas the mean does change per segment.

% Determines GMM components
q(1) = Eval_q(t, phi1);  % the expectation for a full target in segment #1
q(2) = Eval_q(t, phi2);  % the expectation for a full target in segment #2

% Generates samples 
z = GenerateSamples(q, 1);  % a is a function, z(p)

% Builds GMM
fN = @(z) normpdf(z);
cN = @(z) normcdf(z);
fT = @(z,p) nSeg * normpdf(z - p*q');
cT = @(z,p) nSeg * normcdf(z - p*q');

% Builds "th" inverse function, i.e. eta(th) 
icN = @(th) norminv(th);


function g = fevalInRange(f, th, thMin, thMax)
%FEVALINRANGE evaluates g=f(th) for each th, but set INF for off-bounds cases.
g = zeros(size(th));

kTooLow = (th <= thMin);
g(kTooLow) = -inf;

kTooHigh = (th >= thMax);
g(kTooHigh) = inf;

kGood = ~kTooLow & ~kTooHigh;
g(kGood) = arrayfun(f, th(kGood));


function zp = GenerateSamples(q, s)
%GENERATESAMPLES generates samples that well covers both fN and fT.
% The samples are generally [-3*STD : 0.1*STD : MEAN + 3*STD]. Still it considers the case
% of GMM with several means or STDs, and it makes sure a pass through z=0. The return 
% samples are function of the target's power (p), which might affect the MEAN value. 
%
dz = min(s) * 0.1; % fix step as to get 10 samples per sigma, i.e. 60 samples within ±3*sig 
mm1 = -3*max(s);
mm2 = @(p)p*max(q) + 3*max(s);
zp = @(p)[fliplr(0:-dz:mm1), dz:dz:mm2(p)]; % ensure to pass through zero

