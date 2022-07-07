function Show_Figure13a()
%SHOW_FIGURE13A shows Figure 13(a) - Influence of inhomogeneity for different thresholds on synthetic data.

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
numSteps = 101;
maxK = 4;
minK = maxK / (numSteps-1);
ar = 120000; % high as to allow big K values (although 2.01*maxK should be enough)
SNR = 100;  

% sampling ranges
th_all = 10.^(-10:2:-2);
Kb = 10.^([-minK,linspace(minK,maxK,numSteps-1)]); % "-minK" was added for colsing the discontinuity at 0

% layout primitives
t = [SNR,0]'; 
phi0 = SIM_Get_phi(ar);

% init. loop
Ntests = numel(th_all);
Mtests = numSteps;
[pMax, BMax] = deal(zeros(Ntests,Mtests));

% for each "th" calculates pMax(Kb) and Bmax(Kb)
for n = 1:Ntests
    fprintf('test %d/%d\n', n, Ntests);
    th = th_all(n);
    
    % calculate MaxB, MaxP for each K
    for m = 1:numSteps
        Kb0 = Kb(m);
        [~, AGnom] = SIM_Calc_A(phi0, Kb0, 0, t, th, 'Global');
        [~, ALnom] = SIM_Calc_A(phi0, Kb0, 0, t, th, 'Local');
        [phi1, phi2] = SIM_Split_phi(phi0, Kb0, 0);  % extracts the data background
        [~,~,~,~,~,icN] = SIM_Calc_PDF(phi1, phi2, t, 0.5, 'Local');
        etaL = -icN(th);  % inverse function th=>etaL (local eta)
        [~, qL] = SIM_Calc_mu(phi1, phi2, t, 1, 0.5, 'Local');
        p1Exp = log10(min(etaL ./ qL));
        pMax(n,m) = 10.^fminbnd(@(p)-ALnom(10.^p)./AGnom(10.^p), p1Exp-1, p1Exp+1);
        BMax(n,m) = ALnom(pMax(n,m)) ./ AGnom(pMax(n,m));
    end
end
fprintf('done!\n\n');

% init. plots
xAxisRng = [1, max(Kb(:))];  % K([1 end]) for sanity-check
lineWidth = 1.5;
axMrg_ax4 = @(ax)([ax(1:3), 10^(log10(ax(4))+diff(log10(ax(3:4)))*0.02)]);
axMrg = @(ax) axMrg_ax4([ax(1:2), 10^floor(log10(ax(3))), ax(4)]);
vec2legendcell = @(v) split(strtrim(sprintf('th=%.2g ', v)),' ');

% plots Bmax(Kb)
hFig1 = figure; 
loglog(Kb(:), squeeze(BMax(:,:)), 'LineWidth', lineWidth);
axis(axMrg([xAxisRng, 1, max(BMax(:))^1.15]));
xlabel('Kb'); ylabel('Bmax'); grid on;
legend(vec2legendcell(th_all),'location','NorthWest');
set(hFig1,'Position',[680,665,424,313])

% plots pMax(Kb)
hFig2 = figure; 
loglog(Kb(:), squeeze(pMax(:,:)), 'LineWidth', lineWidth);
axis(axMrg([xAxisRng, minmax2(pMax)]));
xlabel('Kb'); ylabel('pmax'); grid on;
legend(vec2legendcell(th_all),'location','NorthEast');
set(hFig2,'Position',[1107,665,424,313]) 
