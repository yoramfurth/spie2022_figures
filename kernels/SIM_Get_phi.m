function phi = SIM_Get_phi(ratio)
%SIM_GET_PHI determines a simple 2D covariance matrix.
%
%Description: 
%    The input of this function refers to a 2D canonic ellipse of r1=1, r2=1/ratio.
%    This determines a 2x2 covariance matrix of [1 0; 0 1/ratio^2]. This simple 
%    covariance is a important cornerstone for building the synthetic data layout,
%    based on linear transforms such as scaling and rotation. 
%    
%Inputs: 
%    ratio - The axes' ratio of a canonic ellipse.
%    
%Outputs:
%    phi - The corresponding covariance matrix.
%
%Example: 
%    phi = SIM_Get_phi(2);  % gets [1 0; 0 0.25]
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

cov11 = 1;
cov22 = cov11 / (ratio^2);
phi = [cov11 0; 0 cov22]; % covariance matrix with no correlation
