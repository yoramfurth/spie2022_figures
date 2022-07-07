function [A, Anom, dAnom] = Eval_A(Pfa, Pd, th)
%EVAL_A evaluates how much relative superior detection is achieved over a range of thresholds.
%
%Description: 
%    This function evaluates a score for the ROC-curve within a range of Pfa values, denoted 0:th (after threshold).  
%    This gives a scalar between 0 and 1, which the higher it is, the higher the ROC-Curve, meaning better performance.
%    The basic idea is to calculate the area gained under the ROC-curve comparing to the null case of a unit function.
%    Its exact expression is: A(th)=(AUC(0:th)-0.5th^2)/(th-0.5th^2), while AUC is the area under the ROC curve. 
%    This implementation has the capability of processing lots of "th" values at once. 
%    More details are available in the paper Section 1.4, and in the thesis report Section 2.6.
% 
%Inputs: 
% 	 Pfa, Pd - Ordered pairs representing the ROC curve.
%    th - The decision thresholds in terms of Pfa. "A" is evaluated within [0:th] range. "th" can be an array.
% 
%Outputs: 
% 	 A - Evaluation of the detection algorithm, A(th), for any of the requested "th".
% 	 Anom - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 
%Example: 
%    th = 0.01;
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    [MF_N, MF_T] = Calc_MF(y, t, 'Global', [], 0.05);  % apply the NGMF
%    [Pfa,Pd] = roc_fast(MF_N, MF_T);  % calculates the ROC curve
%    A = Eval_A(Pfa, Pd, th);
%    fprintf('A(th)=%g for a threshold of th=%g\n', A, th);
% 
%See also ROC_FAST
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
maxIdx = max(2, min([find(Pfa>max(th),1,'first'),numel(Pfa)]));  % finds the maximal Pfa needed to be processed
index_th = (1 : maxIdx); % indexes up to that maximal Pfa
Pfa = Pfa(index_th);
Pd  = Pd(index_th);

% Approximates the Area Under Curve (AUC)
CumA = cumtrapz(Pfa, Pd);  % trapezoid approximation
k = [true, (diff(Pfa)~=0)];  % ignores non-unique points
AUC = interp1(Pfa(k), CumA(k), th, 'linear', 'extrap');  % approximation using linear interpolation

% Evaluates A(th)
Anom = AUC - 0.5*th.^2;  % numerator of A(th)
Adenom = th - 0.5*th.^2;  % denominator of A(th)
Anom = max(0, min(Adenom, Anom));  % fix numerical errors
A = Anom./Adenom;

% Derivative of Anom
if (nargout>=3)	
    Pd_th = interp1(Pfa(k), Pd(k), th, 'next', 'extrap');  % assumes Pfa(end-1)==0, and Pd(end-1)==Pd(end)
    dAnom = Pd_th - th; 
    dAnom = max(0, min(1-th, dAnom));  % fix numerical errors
end
