function Show_Figure4()
%SHOW_FIGURE4 shows Figure 4 - Tests range of effective target's powers through different domains.

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

%% init.
% parameters
p = 0.02;
th = 0.01; 
K = 8;  
theta = 0; 
ar = 120; 
SNR = 100;
badB = 1.1;
badAL = 0.1;

% sampling ranges
pAll = 10.^(-4:0.01:-0.5);

% layout primitives
t = [SNR,0]';  % like snr in SIC data
phi0 = SIM_Get_phi(ar);

% plotting constants
zLim = [-3, 20];
pdfGLim = [0, 0.7];
pdfLLim = [0, 0.42];
pLim = pAll([1,end]);
lineWidth1 = 1.5;
lineWidth2 = 2.5;
lineStyle = '--';
colors0 = get(groot,'defaultAxesColorOrder');
colors = [colors0(1:2,:);       % blue for fN, orange for fT
          colors0(3,:);         % orange for AL<0.1
          [0.42, 0, 0.56];     % purple, for eta,B  
          [0.82, 0.24, 0.52];  % magenta, for Global
		  [0.466 0.674 0.188]]; % green, for Local

%% calculations
% calculcates PDFs at p, and performance (AG,AL,B) at pAll
[AG, AGnom, AGnomLop, fGN, fGT, zG, etaG, muGMax, pGAnc] = SIM_Calc_A_Per_p(pAll, phi0, K, theta, t, p, th, 'Global');
[AL, ALnom, ALnomLop, fLN, fLT, zL, etaL, muLMax, pLAnc] = SIM_Calc_A_Per_p(pAll, phi0, K, theta, t, p, th, 'Local');
B = Eval_B(AG, AL, AGnom, ALnom, AGnomLop, ALnomLop, 1e-14, 0.01);
Bp = interp1(pAll, B, p, 'linear');

% splits B to good and bad regions
kBadB = (B<badB-eps);  % points where B<1
kBadAL = (AL<badAL) & ~kBadB;  % points where AL<0.1
kGoodB = ~(kBadB | kBadAL);  % good points
kGoodB = imfilter(kGoodB, [1 1 1]);  % dilate mask
kBadB = imfilter(kBadB,[1 1 1]) & (kBadB | kBadAL);  % dilate mask
kLines = [kGoodB; kBadB; kBadAL];  % rough distinct lines

%% plots
hFig = figure; set(gcf,'Name','Demo Target-Size Impact'); 
set(hFig,'position',[301 150 1606 836]);

% plots the NGMF PDFs
hAx = subplot(2,2,1);
hp(1) = plot(zG, fGN, 'Color', colors(1,:), 'linewidth', lineWidth1);
grid on; hold on;
hp(2) = plot(zG, fGT, 'Color', colors(2,:), 'linewidth', lineWidth1);
kPd = find(zG>=etaG);
fill(zG(kPd([1,1:end,end])),[0,fGT(kPd),0],colors(2,:),'EdgeColor','none','FaceAlpha',0.3);
plot(etaG, interp1(zG, fGN, etaG, 'linear'), 'o', 'MarkerFaceColor', colors(4,:),'Color', colors(4,:), 'MarkerSize',8);
hp(3) = plotline(etaG, hAx, lineStyle, lineWidth2, colors(4,:));
hp(4) = plotline(muGMax, hAx, lineStyle, lineWidth2, colors(3,:));
axis([zLim, pdfGLim])
xlabel('MF score'); ylabel('probability density'); 
legend(hp, 'without target','with target','\eta_G','p\cdotq_G','location','NorthEast');

% plots the NSMF PDFs
hAx = subplot(2,2,3);
hp(1) = plot(zL, fLN, 'Color', colors(1,:), 'linewidth', lineWidth1);
grid on; hold on;
hp(2) = plot(zL, fLT, 'Color', colors(2,:), 'linewidth', lineWidth1);
kPd = find(zL>=etaL);
fill(zL(kPd([1,1:end,end])),[0,fLT(kPd),0],colors(2,:),'EdgeColor','none','FaceAlpha',0.3);
plot(etaL, interp1(zL, fLN, etaL, 'linear'), 'o', 'MarkerFaceColor', colors(4,:),'Color', colors(4,:), 'MarkerSize',8);
hp(3) = plotline(etaL, hAx, lineStyle, lineWidth2, colors(4,:));
hp(4) = plotline(muLMax, hAx, lineStyle, lineWidth2, colors(3,:));
axis([zLim, pdfLLim])
xlabel('MF score'); ylabel('probability density'); 
legend(hp, 'without target','with target','\eta_L','p\cdotq_1','location','NorthEast');

% plots AG(p) and AL(p)
hAx = subplot(2,2,2);
hP(1) = semilogx(pAll, AG, 'Color', colors(5,:), 'LineWidth', lineWidth1);
grid on; hold on;
hP(2) = semilogx(pAll, AL, 'Color', colors(6,:), 'LineWidth', lineWidth1);
hP(3) = plotline(pGAnc, hAx, lineStyle, lineWidth2, colors(4,:));
plotline(pLAnc, hAx, lineStyle, lineWidth2, colors(4,:));
plot(p, interp1(pAll, AG, p, 'linear'), 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
plot(p, interp1(pAll, AL, p, 'linear'), 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
hP(4) = plotline(p, hAx, lineStyle, lineWidth2, colors(3,:));
axis([pLim, -0.02, 1.02])
xlabel('p'); ylabel('A'); 
legend(hP, 'A_G(p)','A_L(p)','p=\eta/q','p','location','NorthWest')

% plots B(p)
hAx = subplot(2,2,4);
colorsB = colors([4,3,3],:);
lineStyleB = {'-','-','-'};
hP = [];
for k = find(any(kLines,2))'
    toplot = B;
    toplot(~kLines(k,:)) = nan;
    hP(k) = semilogx(pAll, toplot, lineStyleB{k}, 'Color', colorsB(k,:), 'LineWidth', lineWidth1); %#ok
    hold on;
end
BLim = imscale([-0.02 1.02], [0 1], minmax2(B));
hP(3) = plotline(pGAnc, hAx, lineStyle, lineWidth2, colors(4,:));
plotline(pLAnc, hAx, lineStyle, lineWidth2, colors(4,:));
plot(p, Bp, 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
hP(4) = plotline(p, hAx, lineStyle, lineWidth2, colors(3,:));
axis([pLim, BLim])
xlabel('p'); ylabel('B'); grid on;
legend(hP, 'Good B', 'bad B','p=\eta/q','p','location','NorthWest')


function [A, Anom, AnomLop, fN, fT, z, eta, muMax, pAnc] = SIM_Calc_A_Per_p(pAll, phiG, K, theta, t, p, th, MFtype)
%SIM_CALC_A_PER_P calculates A(p) function and related for a set of "p" values.
%
%Description: 
%    Note that this function if a variation of SIM_Calc_A().
%
%Inputs: 
% 	 pAll - A set of "p" values where to evaluate A(p).
% 	 phiG - The global covariance matrix.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 	 t - A vector representing a target's spectrum.
%    p - An exemplar target's power used for creating the auxiary data (PDFs etc.). 
% 	 th - Decision threshold in terms of Pfa.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 A - Evaluation of the detection algorithm, A(th), for any of the requested "p".
% 	 Anom - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 	 fN - PDF of the scores of NSMF with no target.
% 	 fT - PDF of the scores of NSMF with target.
% 	 z - Ticks along the discriminant axis.
% 	 eta - Decision threshold on the scores axis.
%    muMax - The highest expectation in this domain (local or global). 
%    pAnc - The "p" that gives muMax==eta.
% 
%See also SIM_CALC_A
% 

% Splits data
[phi1, phi2] = SIM_Split_phi(phiG, K, theta);  % extracts the data background

% Calculates MF distributions
[fNz, fTzp, zp, cNz, cTzp, icN] = SIM_Calc_PDF(phi1, phi2, t, 0.5, MFtype);
z = zp(p);
fN = fNz(z);
fT = fTzp(z,p);

% Finds special locations (p-anchors)
eta = -icN(th);  % translates th=>eta
[q0, q12] = SIM_Calc_mu(phi1, phi2, t, 1, 0.5, MFtype);
q = max([q0; q12]);
muMax = q*p;
pAnc = eta / q; 

% Calculates A(pAll)
[A, Anom, AnomLop] = SIM_Eval_A(cNz, cTzp, fNz, eta);
A = A(pAll);
Anom = Anom(pAll);
AnomLop = AnomLop(pAll);


