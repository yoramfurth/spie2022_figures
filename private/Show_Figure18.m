function Show_Figure18(data, seg, y, phiG, phiS, E, wvlen)
%SHOW_FIGURE18 shows Figure 18 - The performance of Phi5-eiv1 compared to its most similar pixel.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%    E - Eigenvectors decomposition of phiS{}.
%    wvlen - Wavelengths corresponding to the channels of the data.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
n = 5;
th = 0.01;
sz = size(data);
data2v = reshape(data, [sz(1)*sz(2),sz(3)]); % vectorized data, as is
t1 = E{n}(:,1);  % spectum #1 is the phi5's major eigenvector
corrDataEn = uvec(data2v) * t1;  % correlation with the data
t2 = uvec(data2v(argmax(corrDataEn),:)');  % selects the most correlated pixel

% shows Figure 18(a) - shows the two spectral signatures
figure; set(gcf,'Name','t signature');
tData = [t1,t2];
plot(wvlen, tData);
axis([minmax2(wvlen), minmax2(tData)*1.04]);
grid on; 
xlabel('wavelength'); ylabel('radiance (normalized)')
legend('\Phi5-eiv1','nearest x','Location','NorthEast');
set(gcf,'position',[1360,500,560,280])

% compares statistics
corr12 = max(corrDataEn); 
Kb12 = SIM_Calc_Kb(phiG, phiS, [t1,t2]);
fprintf('corr=%g, ang=%g°, Kb1=%g, Kb2=%g, pos=%d\n', corr12, acosd(corr12), Kb12(1), Kb12(2), all(t2>=0));

Bmax12(1) = Calc_Bmax(y, seg, 1, 0, t1, th, 1);
Bmax12(2) = Calc_Bmax(y, seg, 1, 0, t2, th, 1);
fprintf('Bmax1=%g, Bmax2=%g\n', Bmax12(1), Bmax12(2));

% shows Figure 18(b) - Bmax as function of the angle 
figure; set(gcf,'Name','Bmax(alpha)');
a = [-8:-2, -1.5:0.25:-1, -0.9:0.1:0.9, 1:0.25:1.5, 2:8];
Show_Bmax_XSection(y, seg, t1, t2, a, th, 1, 0, 0);
set(gcf,'position',[1360,0,560,420])
