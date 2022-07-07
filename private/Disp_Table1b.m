function Disp_Table1b(data, seg, y, phiG, phiS, E)
%DISP_TABLE1B displays Table 1(b) - Cooke's cube performance obtained with the optimal target within the positive cone.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%    E - Eigenvectors decomposition of phiS{}.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Finds the best t
avgSNR = Eval_SNR(data, y);  % global SNR
tmax = -SIM_Calc_tmax(phiG, phiS);  % the "minus" performs better (based on offline experiments)

% Defines cross-section axes
n = 4;
eiv1pos = E{n}(:,1)*sign(E{n}(1,1));  % major-eiv, flipped as to have positive 1st component
yax1 = uvec(eiv1pos);  % unit vector in direction of eiv1
yax2 = uvec(tmax);  % unit vector in direction of tmax

% Finds constraint optimum
yax2n = uvec(yax2 - yax1*(yax1'*yax2));  % orthonormal projection
a_rng = Range_Positive_Combs(yax1, yax2n);  % bounds of the constraints area
ta_rng = max(0, sind(a_rng) .* yax2n + cosd(a_rng) .* yax1);  % t at these bounds
Kb_rng = SIM_Calc_Kb(phiG, phiS, ta_rng);  % Kb at these bounds
tmaxCon = ta_rng(:,argmax(Kb_rng));  % the best of them

% Summarizes
disp('tmax at t>0')
tCon = tmaxCon * (avgSNR / Eval_q(tmaxCon, phiG));  % keeps the target's SNR like the global SNR
Disp_Performance_Summary(y, seg, tCon);
