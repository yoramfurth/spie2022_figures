function Show_Figure16(y, seg, phiG, phiS, E)
%SHOW_FIGURE16 shows Figure 16 - segmentation optimum for constraint combinations of tmax with Phi4-eiv1.
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
n = 4;
tMajor = -SIM_Calc_tmax(phiG, phiS);  % the "minus" performed better (based on offline experiments)
eiv1pos = E{n}(:,1)*sign(E{n}(1,1));  % major-eiv, flipped as to have positive 1st component
yax1 = uvec(eiv1pos);  % unit vector in direction of eiv1
yax2 = uvec(tMajor);  % unit vector in direction of tmax

% plots Bmax(a) along the cross-section plane
figure; set(gcf,'Name','Bmax(alpha)');
a = linspace(-180, 180, 256);
[~,tmax] = Show_Bmax_XSection(y, seg, yax1, yax2, a, th, 1, 1, 1);
title(sprintf('\\Phi%d x tmax',n));
set(gcf,'position',[1360,0,560,420])

% shows ellipsoids cross-sections
figure; set(gcf,'Name','Ellipses all');
Show_Ellipses_XSection(phiG, phiS, yax1, yax2, tmax, 1);
xlabel(sprintf('\\Phi%d-eiv%d',n,1))
ylabel('tmax')
set(gcf,'position',[1360,500,560,500])
