function [MF_N, MF_T] = Calc_MF(y, t, MFtype, seg, p)
%CALC_MF calculates the scores of the normalized segmented matched filter (NSMF) without and with target.
%
%Description: 
% 	 The NSMF is a type of matched-filter that deals with segmentation assuming target additive model (see details in readme.txt).
%    This function has two modes: "Global", which calculates the MF globally, and "Local", which calculates it per segment. 
%    Note that the latter requires to provide also a segmentation map. 
% 
%Inputs: 
% 	 y - The residual data to process, generally formed by subtracting the local average from each pixel (x-m).
% 	 t - A vector representing the spectrum of the target of interest. 
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
%    seg - Segmentation map. Required in the "Local" mode only. 
%    p - The power of the target. If empty MF_T outputs a function of "p".
% 
%Outputs: 
% 	 MF_N - The NSMF in the case of no target.
% 	 MF_T - The NSMF with target. If "p" is empty it returns a handle to MF_T(p) function.
% 
%Example: 
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    seg = kmeans_robust(data, 5);  % robustly segments data into 5 clusters
%    [NSMF_N, NSMF_T] = Calc_MF(y, t, 'Local', seg, 0.05);  % apply the NSMF
%    [Pfa,Pd] = roc_fast(NSMF_N, NSMF_T);  % calculates the ROC curve
%    figure; plot(Pfa, Pd);
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% selects between global and local mode
if (strcmpi(MFtype,'Global') || isempty(seg))
    [MF_N, MF_T] = Calc_NGMF(y, t);
else %if MFtype=='Local'
    [MF_N, MF_T] = Calc_NSMF(y, seg, t);
end

if (~isempty(p))
    MF_T = MF_T(p); 
end


function [NGMF_N, NGMF_T] = Calc_NGMF(y, t)
%CALC_NGMF calculates the MF globally

% init.
NGMF_N = zeros(size(y,1),1);  

% evaluates the global matched filter
[NGMF_N(:), NGMF_T_Bias] = Eval_NGMF(y, t);

% applies the bias for the with-target case
NGMF_T = @(p) NGMF_N + p * NGMF_T_Bias; 


function [NSMF_N, NSMF_T] = Calc_NSMF(y, seg, t)
%CALC_NSMF calculates the MF per segment

%init.
seg = seg(:);
sMax = max(seg);
NSMF_N = zeros(size(y,1),1);
NSMF_T_Bias = zeros(sMax,1);

% evaluates the matched filter per segment 
for s=1:sMax
    k = (seg==s);  % mask of the current segment
    yS = y(k,:);  % residual noise of this segment's pixels 
	[NSMF_N(k(:)), NSMF_T_Bias(s)] = Eval_NGMF(yS, t);
end

% applies the bias for the with-target case
NSMF_T = @(p) NSMF_N + p * NSMF_T_Bias(seg(:));  % each pixel gets an addition according to the segment it belongs to


function [MF_N, MF_T_Bias] = Eval_NGMF(y, t)
%EVAL_NGMF evaluates the MF globally

% init.
phi = cov(y);  % estimates covariance from the residual data 
q = Eval_q(t, phi); % the expectation for a full target 

% Applies matched filter without target
MF_STD = q;  % by luck, in the NSMF algorithm, the STD of the matched filter is equal to the target's expectation, so we can reuse it
MF_N = t' * (phi \ y') / MF_STD;  % the NSMF in the case of no target

% The bias the gets the matched filter with target
MF_T_Bias = q; % the bias of the expectation with target in the special case of p==1

