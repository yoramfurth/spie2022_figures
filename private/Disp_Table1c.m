function Disp_Table1c(data, seg, y, phiG, phiS)
%DISP_TABLE1C displays Table 1(c) - Cooke's cube performance obtained with the optimal target that belongs to the data.
%
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map. 
% 	 y - The residual data, generally formed by subtracting the local average from each pixel (x-m).
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Init.
th = 0.01; 

% Calculates Kb at any t=pixel
sz = size(data);
avgSNR = Eval_SNR(data, y);  % global SNR
t = uvec(reshape(data, [sz(1)*sz(2),sz(3)]))'; % the whole data vectorized and normalized
Kb = SIM_Calc_Kb(phiG, phiS, t);

% Finds the pixel whose direction has the highest segmentation benefit
Kb2d = reshape(Kb, sz(1:2));
cc = bwconncomp(Kb2d>2.7); % split to blobs
kBigKb = cellfun(@(idx)idx(argmax(Kb(idx))), cc.PixelIdxList); % the highest in each blob (just time saving)
Ntests = numel(kBigKb);
BMax = deal(zeros(Ntests,1));
for n=1:Ntests
    k = kBigKb(n);
    [yMaxKb, xMaxKb] = ind2sub(sz(1:2),k);
    tData = squeeze(data(yMaxKb,xMaxKb,:));  % said to be max, let's see
    t = tData * (avgSNR / Eval_q(tData, phiG));  % keeps the target's SNR like the global SNR

    [~, qL, th2etaL] = Calc_mu(y, seg, 1, 0, t, 1, 'Local');
    p1 = min(th2etaL(th) ./ qL);  % this is "p1" anchor which is the smallest local bound
    BMax(n) = Calc_B(y, seg, 1, 0, t, p1, th);  % the real Bmax is close to "p1" anchor
end
nBest = argmax(BMax);
[yMaxKb, xMaxKb] = ind2sub(sz(1:2),kBigKb(nBest));

% Summarizes
disp('tmax in data')
tmax = squeeze(data(yMaxKb,xMaxKb,:));
t = tmax * (avgSNR / Eval_q(tmax, phiG));  % keeps the target's SNR like the global SNR
Disp_Performance_Summary(y, seg, t);
