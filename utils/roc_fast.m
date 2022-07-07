function [Pfa,Pd] = roc_fast(F1,F2)
%ROC_FAST Computes the Receiver Operating Characteristic (ROC) curve
%
%Description: 
%    This function computes the ROC curve given two 1D decision sets, 
%    e.g. along a discriminant axis. One set is related to the data 
%    without target, from which the false-alarms (Pfa) are extracted. 
%    The second set is related to the data with target, which defines 
%    the probability of detection (Pd) as function of the respective 
%    Pfa. The resulted curve is determined by pairs of (Pfa,Pd), which 
%    are outputted as 2 separated vectors.  
%
%    Note that this function is equivalent to perfcurve() but performs 
%    4 times faster. 
%    
%Inputs: 
% 	 F1 - A decision set related to data without target. 
%    F2 - A decision set related to data with target. 
% 
%Outputs: 
%    Pfa - Probability of false alarm
%    Pd - Probability for detection
%
%Example: 
%    F1 = SomeMatchedFilter(x);
%    F2 = SomeMatchedFilter(x + t);
%    [Pfa, Pd] = roc_fast(F1, F2);
%    figure; plot(Pfa, Pd);
% 
%See also: ROC, PERFCURVE
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% prepare LUTs
F1 = F1(:)';
F2 = F2(:)';
s1 = sort(F1,'descend');
s2 = sort(F2,'descend');
s12 = sort([F1,F2],'descend');

t1 = (1:numel(s1))/numel(s1);
t2 = (1:numel(s2))/numel(s2);

% extrapolation
if (s1(1)<s2(1))
    s1 = [s2(1), s1];
    t1 = [0, t1];
elseif (s1(1)>s2(1))
    s2 = [s1(1), s2];
    t2 = [0, t2];
end
if (s1(end)<s2(end))
    s2 = [s2, s1(end)];
    t2 = [t2, 1];
elseif (s1(end)>s2(end))
    s1 = [s1, s2(end)];
    t1 = [t1, 1];
end

% calculate ROC
Pd = [0, interp1_unique(s2,t2,s12,'next')];
Pfa = [0, interp1_unique(s1,t1,s12,'next')];


function Vout = interp1_unique(X,V,varargin)
% Bypasses "The grid vectors must contain unique points" error
[~, ind] = unique(X);  % ind = index of first occurrence of a repeated value
Vout = interp1(X(ind), V(ind), varargin{:});
