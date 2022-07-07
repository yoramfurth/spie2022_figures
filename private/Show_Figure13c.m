function Show_Figure13c(data, seg, y, phiG, phiS, E)
%SHOW_FIGURE13C shows Figure 13(c) - Influence of inhomogeneity for different thresholds on Cooke's data-cube.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%    E - Eigenvectors decomposition of phiS{}.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
Mtests = 40;  % number of eigenvectors to proceed in each segment
th = 0.01;
avgSNR = Eval_SNR(data, y);  % global SNR

% init. loop
Ntests = numel(phiS);
[Kb, BMax, pMax] = deal(zeros(Ntests, Mtests));

% for each eigenvector calculates Kb, pMax, Bmax
for n = 1:Ntests
    Tn = E{n};
    for m = 1:Mtests
        fprintf('Ntest %d/%d => Mtest %d/%d\n', n, Ntests, m, Mtests);  % takes time!
        t = Tn(:,m);
        Kb(n,m) = SIM_Calc_Kb(phiG, phiS, t);
        t = t * (avgSNR / Eval_q(t, phiG));  % keeps the target's SNR like the global SNR
        [BMax(n,m), pMax(n,m)] = Calc_Bmax(y, seg, 1, 0, t, th);
    end
end

% plots a scattergram of Bmax(Kb)
xAxisRng = [1, max(Kb(:))];  
hFig1 = figure; 
loglog(Kb(:), BMax(:), '.');
axis; axis([xAxisRng, 1, max(BMax(:)).^1.15])
xlabel('Kb'); ylabel('Bmax'); grid on;
legend('th=0.01','location','SouthEast');
set(hFig1,'Position',[680,665,424,313])

% plots a scattergram of pMax(Kb)
hFig2 = figure; 
loglog(Kb(:), pMax(:), '.');
axis([xAxisRng, 0.05, max(pMax(:))])
xlabel('Kb'); ylabel('pmax'); grid on;
legend('th=0.01','location','NorthEast');
set(hFig2,'Position',[1107,665,424,313])  

