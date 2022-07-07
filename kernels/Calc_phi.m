function [phiG, phiS] = Calc_phi(y, seg)
%CALC_PHI estimates the local and global covariances of a given data.
%
%Description: 
% 	 A covariance calculates the correlation of any pair of "yi-E(yi)", while E() is
%    expectation (estimated by averaging), and "yi" is one column, i.e. y(:,i). Note 
%    that since "y" is normally after subtracting the average, and internal E(yi) gives
%    about zero, so the operation is almost equivalent to simply y'*y/N. The local 
%    covariances refers to a calculation per segment, defined by "seg" segmentation map.
% 
%Inputs: 
% 	 y - Data to process, normally x-m.
%    seg - Segmentation map.
% 
%Outputs: 
% 	 phiG - Global covariance matrix.
% 	 phiS - Cell array of local covariance matrixes.
% 
%See also COV
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

phiG = cov(y,1);  % "1" because this gives an estimation that better matches our case
phiSegFun = @(s) cov(y(seg(:)==s,:),1);  
phiS = arrayfun(phiSegFun, unique(seg(:)), 'UniformOutput', false);

