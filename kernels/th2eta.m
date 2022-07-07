function eta = th2eta(MF_N, th)
%TH2ETA Converts the decision threshold from the Pfa domain to the scores domain.
%
%Description: 
%    This function calculates the inverse tranform from the Pfa domain (given by cumulative
%    false alarms, and denoted "th") back to the scores domain.
% 
%Inputs: 
% 	 MF_N - The NSMF in the case of no target.
% 	 th - Decision threshold in terms of Pfa.
% 
%Outputs: 
% 	 eta - Decision threshold on the scores axis.
% 
%Example: 
%    th = 0.01;
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    MF_N = Calc_MF(y, t, 'Global', [], 0.05);  % apply the NSMF
%    eta = th2eta(MF_N, th);  % calculates the ROC curve
%    fprintf('eta=%g is the threshold to gets Pfa=%g\n', eta, th);
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

eta = prctile(MF_N, 100*(1-th));
