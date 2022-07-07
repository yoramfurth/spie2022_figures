function Disp_Table1a(data, seg, y, phiG, phiS)
%DISP_TABLE1A displays Table 1(a) - Cooke's cube performance obtained with the optimal target under no constraints.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Finds the best t
avgSNR = Eval_SNR(data, y);  % global SNR
tmax = -SIM_Calc_tmax(phiG, phiS);  % the "minus" performs better (based on offline experiments)

% Summarizes
fprintf('tmax')
t = tmax * (avgSNR / Eval_q(tmax, phiG));  % keeps the target's SNR like the global SNR
Disp_Performance_Summary(y, seg, t); 
