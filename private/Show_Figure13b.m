function Show_Figure13b(SIC_data, SIC_seg, t)
%SHOW_FIGURE13B shows Figure 13(b) - Influence of inhomogeneity for different thresholds on SIC data-cube.
%
%Inputs: 
% 	 SIC_data - SIC data-cube. 
% 	 SIC_seg - SIC data-cube's segmentation (2 segments). 
%    t - A vector that represents a target spectrum.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
maxK = 4;
numSteps = 101;
minK = maxK / (numSteps-1);
SNR = 100; 

% sampling ranges
th_all = 10.^(-3:0.5:-1);
Kb = 10.^([-minK,linspace(minK,maxK,numSteps-1)]); % "-minK" was added for closing a discontinuity at 0
assert(all(Kb~=1),'Kb==1 is not supported');

% init. loop
Ntests = numel(th_all);
Mtests = numSteps;
[pMax, BMax] = deal(zeros(Ntests,Mtests));

% for each "th" calculates pMax(Kb) and Bmax(Kb)
for m = 1:Mtests
    fprintf('test %d/%d\n', m, Mtests);
    K = Kb(m);  % scaling factor
    actualSNR = Calc_mu(SIC_data, SIC_seg, K, 0, t, 1, 'Global');
    tNorm = t * (SNR / actualSNR); % defines a target with the requested SNR
    
    % calculates pMax(Kb) and Bmax(Kb)
    for n = 1:Ntests
        th = th_all(n);
        [BMax(n,m), pMax(n,m)] = Calc_Bmax(SIC_data, SIC_seg, K, 0, tNorm, th);
    end
end
fprintf('done!\n\n');

% init. plots
xAxisRng = [1, max(Kb(:))];  % K([1 end]) for sanity-check
lineWidth = 1.5;
axMrg_ax4 = @(ax)([ax(1:3), 10^(log10(ax(4))+diff(log10(ax(3:4)))*0.02)]);
axMrg = @(ax) axMrg_ax4([ax(1:2), 10^floor(log10(ax(3))), ax(4)]);
vec2legendcell = @(v) split(strtrim(sprintf('th=%.2g ', v)),' ');
m2v = @(m)m(:);
yAxisRngP = minmax2([m2v(pMax(:,:))',4.1e-2,1e-6]); 

% plots Bmax(Kb)
hFig1 = figure; 
loglog(Kb(:), squeeze(BMax(:,:)), 'LineWidth', lineWidth);
axis(axMrg([xAxisRng, 1, max(BMax(:))^1.15]));
xlabel('Kb'); ylabel('Bmax'); grid on;
legend(vec2legendcell(th_all),'location','SouthEast');
set(hFig1,'Position',[680,665,424,313])

% plots pMax(Kb)
hFig2 = figure; 
loglog(Kb(:), squeeze(pMax(:,:)), 'LineWidth', lineWidth);
axis(axMrg([xAxisRng, yAxisRngP]));
xlabel('Kb'); ylabel('pmax'); grid on;
legend(vec2legendcell(th_all),'location','NorthEast');
set(hFig2,'Position',[1107,665,424,313]) 
