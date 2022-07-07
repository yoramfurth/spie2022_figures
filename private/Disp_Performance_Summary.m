function Disp_Performance_Summary(y, seg, t)
%DISP_PERFORMANCE_SUMMARY calculates principal performance metrics and displays a summary.
%
%Inputs: 
% 	 y - The residual data to process, generally formed by subtracting the local average from each pixel (x-m).
%    seg - Segmentation map. 
% 	 t - A vector representing the spectrum of the target of interest. 
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Init.
th = [0.001, 0.01, 0.1];

% Calculates Kb 
[qG, ~, th2etaG] = Calc_mu(y, seg, 1, 0, t, 1, 'Global');
[~, qL, th2etaL] = Calc_mu(y, seg, 1, 0, t, 1, 'Local');
Kb = (max(qL)/qG) .* ones(size(th));

% Calculates dp
pG = th2etaG(th) ./ qG;
p1 = th2etaL(th) ./ max(qL);
dp = pG ./ p1;

% Calculates p-max
pMax = arrayfun(@(th0) Calc_pmax(y, seg, 1, 0, t, th0), th);

% Calculates performance metrics
[B,AG,AL] = Calc_B(y, seg, 1, 0, t, pMax, num2cell(th));

% Displays a summary
TT = array2table([Kb; dp; pMax; AG; AL; B],...
    'VariableNames',{'A001' 'A01' 'A1'},...
    'RowName',{'K_b', 'dp', 'p_max', 'NGMF (A_G)', 'NSMF (A_L)', 'Benefit (B)'});
disp(TT)
