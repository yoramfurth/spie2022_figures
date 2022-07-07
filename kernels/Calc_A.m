function [A, Anom, dAnom] = Calc_A(data, seg, K, theta, t, p, th, MFtype)
%CALC_A calculates A as function of various factors.
%
%Description: 
% 	 The metric "A" evaluates how much relative superior detection is achieved over a range of thresholds.
%    CALC_A receives basic spectral factors and calculates "A" from zero, by computing all modules up to
%    the required point. This function can test set of values for any relevant factor. For example, if "p" 
%    is an array, then the result is "A(p)" for any "p" in the array. 
% 
%Inputs: 
% 	 data (or y) - The data-cube to process (of the residual noise).
%    seg - Segmentation map.
%    K - Total scaling of the two segments. Null is 1. Supports array.
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. Supports array.
% 	 t - A vector representing the spectrum of the target of interest. Supports array.
%    p - The target's power. Supports array.
%    th - The decision thresholds in terms of Pfa. Supports array.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
%
%Outputs: 
% 	 A - Evaluation of the detection algorithm, A(th), for any of the requested "th".
% 	 Anom - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    th = 1./10.^1:0.1:3;  % th from 0.001 to 0.1
%    A = Calc_A(data, seg, 4, 20, t, 1, th, 'Local');  % test spiting into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    figure; semilogx(th, A);  % plots A(th) from 0.001 to 0.1
% 
%See also EVAL_A
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init. variables
if (~iscell(t))
    t = {t};
end
if (~iscell(th))
    th = {th(:)};
end
if (isempty(K))
    K = 1;
end
if (isempty(theta))
    theta = 0;
end

% init. for-loop flowMng
flowMng = FlowManager();  % a tricky class for efficient runtime
flowMng = flowMng.AddLabel(K,theta);
flowMng = flowMng.AddLabel(t);
flowMng = flowMng.AddLabel(p);

% prepare input arrays for loop
mainArg = GetMainArg(K, theta, t, p);
sz = size(mainArg);
K = K.*ones(sz);  % supports vector or scalar
theta = theta.*ones(sz);
t(1:sz(1),1:sz(2)) = t;
p = p.*ones(sz);  % supports vector or scalar

th(1:sz(1),1:sz(2)) = th;
th = [th{:}];
szTh = size(th);

% prepare output arrays
A = zeros(szTh);
Anom = zeros(szTh);
dAnom = zeros(szTh);
NTests = prod(sz);

% calculates performance
y = Calc_y(data, seg, 1, 1);
for nK=1:NTests
    [A(:,nK), Anom(:,nK), dAnom(:,nK)] = ...
        Calc_A_single_test(y, seg, K(nK), theta(nK), t{nK}, p(nK), th(:,nK));
end


    function [A, Anom, dAnom] = Calc_A_single_test(y, seg, K, theta, t, p, th)
        %Label1: (split data)
        flowMng = flowMng.SetLabel(1);
        if (flowMng.start())
            y = Split_y(y, seg, K, theta);  % splits into 2 distinct segments
            flowMng = flowMng.SetReg(y);
        else
            y = flowMng.GetReg();
        end
        
        %Label2: (MF detection)
        flowMng = flowMng.SetLabel(2);
        if (flowMng.start())
            [MF_N, MF_Tp] = Calc_MF(y, t, MFtype, seg, []);
            flowMng = flowMng.SetReg(MF_N, MF_Tp);
        else
            [MF_N, MF_Tp] = flowMng.GetReg();
        end
        
        %Label3: (ROC curve)
        flowMng = flowMng.SetLabel(3);
        if (flowMng.start()) % save TPT while looping on "p"
            MF_T = MF_Tp(p);  
            [Pfa,Pd] = roc_fast(MF_N,MF_T);
            flowMng = flowMng.SetReg(Pfa,Pd);
        else
            [Pfa,Pd] = flowMng.GetReg();
        end           
        
        %Label4: (Calculates Ath)
        [A,Anom,dAnom] = Eval_A(Pfa,Pd,th);
    end
end


function mainArg = GetMainArg(varargin) 
% Get the main argument among the others.
% A "main argument" is an array on which the funciton is evaluated. 
[~, kMax] = max(cellfun(@numel,varargin));  % supports vector
mainArg = varargin{kMax};
end

