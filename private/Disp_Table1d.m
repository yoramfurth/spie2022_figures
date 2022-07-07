function Disp_Table1d(data, seg, t)
%DISP_TABLE1D displays Table 1(d) - Cooke's cube performance obtained with the optimal target that belongs to a given set of real targets.
%
%Inputs: 
% 	 data - A data-cube, given in reflectance units. 
%    seg - Segmentation map. 
%    t - A target's spectrum, sampled by the data authors.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Calculates the residual noise (x-m)
y = Calc_y(data, seg, 1, 0);

% Summarizes
disp('tmax in targets')
Disp_Performance_Summary(y, seg, t);
