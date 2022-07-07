function Show_Figure1(data, seg, t)
%SHOW_FIGURE1 shows Figure 1 - Successful versus failed segmentation on Cooke's data.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 t - A vector representing the spectrum of the target of interest. 
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init
pickRow = @(V,k)V(:,k);
th = 0.01;
y = Calc_y(data, [], 0, 0);  % simple x-m
[~, phiS0] = Calc_phi(y, seg);  % covariances with simple mean

% A good case
tNorm = pickRow(eig_sorted(phiS0{4}),3) * norm(t);
p = 0.015349;  % calculated as to give Ath=0.12 on this "t"
Show_NSMF_BasicGraphs(data, seg, tNorm, p, th); 

% A bad case
tNorm = pickRow(eig_sorted(phiS0{1}),1) * norm(t);
p = 0.052764;  % calculated as to give Ath=0.12 on this "t" 
Show_NSMF_BasicGraphs(data, seg, tNorm, p, th); 


function Show_NSMF_BasicGraphs(data, seg, t, p, th)
%SHOW_NSMF_BASICGRAPHS compares NGMF vs NSMF distributions
%
%Description:
%    This function compares NGMF vs NSMF distributions, given data and target.
%    For each MF it first plots 2 distributions, without target and with target,
%    and then is plots the respective ROC curves.
%
%Inputs: 
% 	 data - The data-cube to process.
%    seg - Segmentation map.
% 	 t - A vector representing the spectrum of the target of interest.
%    p - The target's power.
%    th - The decision thresholds in terms of Pfa.
%

%% Plots PDF of NGMF 
y = Calc_y(data, [], 0, 0);  % simple x-m
[NGMF_N, NGMF_T] = Calc_MF(y, t, 'Global', seg, p);
[hFig1, histMax] = plotPDF(NGMF_N, NGMF_T, th, 'Scores on Global Data');
ax1 = [-6, 18, 0, 1.1*histMax]; axis(ax1);
set(hFig1,'position',[373, 388, 560, 295]);

%% Plots PDF of NSMF
[NSMF_N, NSMF_T] = Calc_MF(y, t, 'Local', seg, p);
hFig2 = plotPDF(NSMF_N, NSMF_T, th, 'Scores on Segmented Data');
axis(ax1);
set(hFig2,'position',[373, 45, 560, 295]);

%% Plots ROC curves
[PfaG, PdG] = roc_fast(NGMF_N, NGMF_T);
[PfaS, PdS] = roc_fast(NSMF_N, NSMF_T);
hFig3 = plotRocCurves(PfaG, PdG, PfaS, PdS);
axis([0 th 0 1]);
set(hFig3, 'position', [934, 224, 601, 460]);


function [hFig, histMax] = plotPDF(MF_N, MF_T, th, ttl)
%PLOTPDF plots 2 distributions, without target and with target.
%
%Description: 
% 	 MF_N - The NSMF in the case of no target.
% 	 MF_T - The NSMF with target. 
%    th - The decision thresholds in terms of Pfa. 
%    ttl - Title of the figure. 
%
%Outputs: 
% 	 hFig - Handle to the displayed figure.
% 	 histMax - The maximum of both distributions.
%
nBins = 4000;
lwid = 1;
th_color = [0.42, 0, 0.56];
ttl_fntSize = 18;

NormHist=@(h,x)h/(sum(h)*mean(diff(x)));
[h1,x1] = hist(MF_N,nBins); h1=NormHist(h1,x1);
[h2,x2] = hist(MF_T,nBins); h2=NormHist(h2,x2);
zTh = prctile(MF_N,100*(1-th));

hFig = figure; set(gcf,'Name',ttl);
plot(x1,h1,x2,h2,'linewidth',lwid); grid on; hold on;
plotline(zTh,gca(hFig),':',3,th_color);
title(ttl,'FontSize',ttl_fntSize)
xlabel('MF score'); ylabel('probability density');
legend('without target','with target','threshold');

histMax = max(max(h1(:)),max(h2(:)));


function hFig = plotRocCurves(PfaG, PdG, PfaS, PdS)
%PLOTROCCURVES plots the global and the local ROC curves.
%
%Description: 
%    PfaG - Probability of NGMF false alarm
%    PdG - Probability for NGMF detection
%    PfaS - Probability of NSMF false alarm
%    PdS - Probability for NSMF detection
%
%Outputs: 
% 	 hFig - Handle to the displayed figure.
%
lwid = 3;
clrG = [0.82, 0.24, 0.52];
clrS = [0.466, 0.674, 0.188];
ttl_fntSize = 18;
lbl_fntSize = 14;
legand_fntSize = 14;

hFig = figure; set(gcf,'Name','ROC curves on Segmented Data');
plot(PfaG, PdG, 'LineWidth', lwid, 'color', clrG); grid on; hold on; 
plot(PfaS, PdS, 'LineWidth', lwid, 'color', clrS)

title('ROC Curves', 'FontSize', ttl_fntSize);
xlabel('P_{FA}', 'FontSize', lbl_fntSize);
ylabel('P_{d}', 'FontSize', lbl_fntSize);
hLeg = legend('Global MF', 'Normalized SMF');
set(hLeg,'FontSize', legand_fntSize, 'Location', 'northwest');
