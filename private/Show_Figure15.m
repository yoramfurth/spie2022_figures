function Show_Figure15(y, seg, phiG, phiS, E)
%SHOW_FIGURE15 shows Figure 15 - segmentation performance for linear combinations of tmax with Phi5-eiv1.
%
%Inputs: 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    seg - Segmentation map. 
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%    E - Eigenvectors decomposition of phiS{}.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% defines axes of a cross-section plane
th = 0.01; 
n = 5;
tMajor = -SIM_Calc_tmax(phiG, phiS);  % the "minus" performs better (based on offline experiments)
yax1 = uvec(tMajor);  % unit vector in direction of tmax
yax2 = uvec(E{n}(:,1));  % unit vector in direction of eiv1

% plots Bmax(a) along the cross-section plane
figure; set(gcf,'Name','Bmax(alpha)');
a = linspace(-180, 180, 64);
[~,tmax] = Show_Bmax_XSection(y, seg, yax1, yax2, a, th, 0, 1, 1);
title(sprintf('tmax x \\Phi%d',n));
set(gcf,'position',[1360,0,560,420])

% shows ellipsoids cross-sections
figure; set(gcf,'Name','Ellipses all');
Show_Ellipses_XSection(phiG, phiS, yax1, yax2, tmax, 0);
xlabel('tmax')
ylabel(sprintf('\\Phi%d-eiv%d',n,1))
set(gcf,'position',[1360,500,560,500])
