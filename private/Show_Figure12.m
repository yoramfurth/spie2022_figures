function Show_Figure12()
%SHOW_FIGURE12 shows Figure 12 - Tests the influence of the error threshold in different domains.

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

%% init.
% parameters
K = 8;  
th = 10^(-2.4);
theta = 0; 
ar = 120; 
SNR = 100;
badB = 1.1;
badAL = 0.1;

% sampling ranges
thAll = 10.^(-8:0.1:-1);

% layout primitives
t = [SNR,0]';  % like snr in SIC data
phi0 = SIM_Get_phi(ar);

% plotting constants
zLim = [-3, 10];
pdfGLim = [0, 0.7];
pdfLLim = [0, 0.42];
lineWidth1 = 1.5;
lineWidth2 = 2.5;
lineStyle = '--';
colors0 = get(groot,'defaultAxesColorOrder');
colors = [colors0(1:2,:);       % blue for fN, orange for fT
          colors0(3,:);         % orange for AL<0.1
          [0.42, 0, 0.56];      % purple, for eta,B  
          [0.82, 0.24, 0.52];   % magenta, for Global
		  [0.466 0.674 0.188]]; % green, for Local

%% calculations
% calculcates PDFs at th, and performance (AG,AL,B) at thAll
[AL, ALnom, ALnomLop, fLN, fLT, zL, etaL, muLMax, p1, etaAll, pAll] = SIM_Calc_A_Per_th(thAll, phi0, K, theta, t, [], th, [], 'Local');
[AG, AGnom, AGnomLop, fGN, fGT, zG, etaG, muGMax] = SIM_Calc_A_Per_th(thAll, phi0, K, theta, t, p1, th, pAll, 'Global');
B = Eval_B(AG, AL, AGnom, ALnom, AGnomLop, ALnomLop, 1e-14, 0.01);
BEta = interp1(etaAll, B, etaL, 'linear');

etaLim = minmax2(etaAll);

% splits B to good and bad regions
kBadB = (B<badB-eps);  % points where B<1
kBadAL = (AL<badAL) & ~kBadB;  % points where AL<0.1
kGoodB = ~(kBadB | kBadAL);  % good points
kGoodB = imfilter(kGoodB, [1 1 1]);  % dilate mask
kBadB = imfilter(kBadB,[1 1 1]) & (kBadB | kBadAL);  % dilate mask
kLines = [kGoodB; kBadB; kBadAL];  % rough distinct lines

%% plots
hFig = figure; set(gcf,'Name','Demo Threshold Impact'); 
set(hFig,'position',[301 150 1606 836]);

% plots the NGMF PDFs
hAx = subplot(2,2,1);
hP = [];
hP(1) = plot(zG, fGN, 'Color', colors(1,:), 'linewidth', lineWidth1);
grid on; hold on;
hP(2) = plot(zG, fGT, 'Color', colors(2,:), 'linewidth', lineWidth1);
kPd = find(zG>=etaG);
if (~isempty(kPd))
    fill(zG(kPd([1,1:end,end])),[0,fGT(kPd),0],colors(2,:),'EdgeColor','none','FaceAlpha',0.3);
end
plot(etaG, interp1(zG, fGN, etaG, 'linear'), 'o', 'MarkerFaceColor', colors(4,:),'Color', colors(4,:), 'MarkerSize',8);
hP(3) = plotline(etaG, hAx, lineStyle, lineWidth2, colors(4,:));
hP(4) = plotline(muGMax, hAx, lineStyle, lineWidth2, colors(3,:));
axis([zLim, pdfGLim])
xlabel('MF score'); ylabel('probability density'); 
legend(hP, 'without target','with target','\eta_G','\mu_G','location','NorthEast');

% plots the NSMF PDFs
hAx = subplot(2,2,3);
hP = [];
hP(1) = plot(zL, fLN, 'Color', colors(1,:), 'linewidth', lineWidth1);
grid on; hold on;
hP(2) = plot(zL, fLT, 'Color', colors(2,:), 'linewidth', lineWidth1);
kPd = find(zL>=etaL);
fill(zL(kPd([1,1:end,end])),[0,fLT(kPd),0],colors(2,:),'EdgeColor','none','FaceAlpha',0.3);
plot(etaL, interp1(zL, fLN, etaL, 'linear'), 'o', 'MarkerFaceColor', colors(4,:),'Color', colors(4,:), 'MarkerSize',8);
hP(4) = plotline(muLMax, hAx, lineStyle, lineWidth1, colors(3,:));
hP(3) = plotline(etaL, hAx, lineStyle, lineWidth2, colors(4,:));
axis([zLim, pdfLLim])
xlabel('MF score'); ylabel('probability density'); 
legend(hP, 'without target','with target','\eta_L','\mu_{max}','location','NorthEast');

% plots AG(eta) and AL(eta)
hAx = subplot(2,2,2);
hP = [];
hP(1) = semilogy(etaAll, AG, 'Color', colors(5,:), 'LineWidth', lineWidth1);
grid on; hold on;
hP(2) = semilogy(etaAll, AL, 'Color', colors(6,:), 'LineWidth', lineWidth1);
plot(etaL, interp1(etaAll, AG, etaL, 'linear'), 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
plot(etaL, interp1(etaAll, AL, etaL, 'linear'), 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
hP(3) = plotline(etaL, hAx, lineStyle, lineWidth2, colors(3,:));
axis([etaLim, min(AG(:))*0.8, 1.02])
xlabel('\eta'); ylabel('A'); 
legend(hP, 'A_G(\eta)','A_L(\eta)','\eta','location','SouthWest')

% plots B(eta)
hAx = subplot(2,2,4);
colorsB = colors([4,3,3],:);
lineStyleB = {'-','-','-'};
hP = [];
linesTypes = {'$$\hat{B}_{max}$$', '$$bad B$$', '$$bad AL$$'};
actualLineTypes = {};
for k = find(any(kLines,2))'
    toplot = B;
    toplot(~kLines(k,:)) = nan;
    hP(end+1) = semilogy(etaAll, toplot, lineStyleB{k}, 'Color', colorsB(k,:), 'LineWidth', lineWidth1); %#ok
    actualLineTypes{end+1} = linesTypes{k};  %#ok
    hold on;
end
BLim = [min(B(:))*0.8, max(B(:))*1.2];
plot(etaL, BEta, 'o', 'MarkerFaceColor', colors(3,:),'Color', colors(3,:), 'MarkerSize',8);
hP(end+1) = plotline(etaL, hAx, lineStyle, lineWidth2, colors(3,:));
axis([etaLim, BLim])
xlabel('\eta'); ylabel('B'); grid on;
legend(hP, actualLineTypes{:},'$$\eta$$','Interpreter','Latex','location','NorthWest')


function [A, Anom, AnomLop, fN, fT, z, eta, muMax, p, etaAll, pAll] = SIM_Calc_A_Per_th(thAll, phiG, K, theta, t, p, th, pAll, MFtype)
%SIM_CALC_A_PER_TH calculates A(p) function and related for a set of "th" values.
%
%Description: 
%    Note that this function if a variation of SIM_Calc_A().
%
%Inputs: 
% 	 thAll - A set of "th" values where to evaluate A(th).
% 	 phiG - The global covariance matrix.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 	 t - A vector representing a target's spectrum.
%    p - An exemplar target's power used for creating the auxiary data. Supports empty on "Local" MFtype.
% 	 th - Decision threshold in terms of Pfa.
% 	 pAll - A set of "p" values where to evaluate A(th,p). Supports empty on "Local" MFtype.
% 	 MFtype - "Global" for NGMF, and "Local" for NSMF.
% 
%Outputs: 
% 	 A - Evaluation of the detection algorithm, A(th), for any of the requested "th".
% 	 Anom - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 	 fN - PDF of the scores of NSMF with no target.
% 	 fT - PDF of the scores of NSMF with target.
% 	 z - Ticks along the discriminant axis.
% 	 eta - Decision threshold on the scores axis that corresponds to "th".
%    muMax - The highest expectation in this domain (local or global). 
%    p - On empty input "p", returns an approximation to pMax at "th". Available on "Local" MFtype only.
%    etaAll - Decision thresholds that correspond to "thAll".
%    pAll - On empty input "pAll", returns an approximation to pMax at "thAll". Available on "Local" MFtype only.
% 
%See also SIM_CALC_A
% 

% init.
pMaxFactor = 1.1;  %sqrt(2);  % factor that approximates "pMax" from "p1"
assert(strcmpi(MFtype,'Local') || ~isempty(p), 'empty p is relevant only in Local mode')
assert(strcmpi(MFtype,'Local') || ~isempty(pAll), 'empty pAll is relevant only in Local mode')

% Splits data
[phi1, phi2] = SIM_Split_phi(phiG, K, theta);  % extracts the data background

% Calculates MF distributions (analytic)
[fNz, fTzp, zp, cNz, cTzp, icN] = SIM_Calc_PDF(phi1, phi2, t, 0.5, MFtype);

% Finds special locations (p-anchors)
eta = -icN(th);  % translates th=>eta
[q0, qL] = SIM_Calc_mu(phi1, phi2, t, 1, 0.5, MFtype);
q = max([q0; qL]);

if (strcmpi(MFtype,'Local') && isempty(p))  % on empty "p" uses pMax approximation
	p1 = eta / q;  % "p1" anchor for this "eta"
	p = p1 * pMaxFactor;  % approximates pMax
end

muMax = q*p;

% Samples MF distributions
z = zp(p);
fN = fNz(z);
fT = arrayfun(@(z0)fTzp(z0,p), z);

% Calculates A(thAll)
etaAll = -icN(thAll);  % translates th=>eta
[A, Anom, AnomLop] = SIM_Eval_A(cNz, cTzp, fNz, etaAll);

if (strcmpi(MFtype,'Local') && isempty(pAll))  % on empty "pAll" uses pMax approximation 
	p1 = etaAll / q; % "p1" anchor for every "eta"
    pAll = p1 * pMaxFactor;  % approximates pMax for every "th" 
else
    assert(numel(pAll)==numel(thAll),'array size should be equal');
end

A = A(pAll);
Anom = Anom(pAll);
AnomLop = AnomLop(pAll);
