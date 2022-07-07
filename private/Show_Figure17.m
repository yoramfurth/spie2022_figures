function Show_Figure17(data, seg, y, phiG, phiS)
%SHOW_FIGURE17 shows Figure 17 - several views with the best performing pixel.
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

% init.
th = 0.01; 
sz = size(data);
avgSNR = Eval_SNR(data, y);  % global SNR

% calculates Kb for each t=pixel
t = uvec(reshape(data, [sz(1)*sz(2),sz(3)]))'; % the whole data vectorized and normalized
Kb = SIM_Calc_Kb(phiG, phiS, t);

% finds the most beneficial pixel as a target
Kb2d = reshape(Kb, sz(1:2));
cc = bwconncomp(Kb2d>2.7); % splits to blobs
kBigKb = cellfun(@(idx)idx(argmax(Kb(idx))), cc.PixelIdxList); % finds the highest Kb in each blob (just time saving)
Ntests = numel(kBigKb);
BMax = deal(zeros(Ntests,1));  % Bmax per pixel
for n=1:Ntests  % for each blob, calculates Bmax
    k = kBigKb(n);  % index of of the highest Kb in this blob
    [yMaxKb, xMaxKb] = ind2sub(sz(1:2),k);  % index to coordinates
    tData = squeeze(data(yMaxKb,xMaxKb,:));  % takes this pixel
    t = tData * (avgSNR / Eval_q(tData, phiG));  % keeps the target's SNR like the global SNR

    [~, qL, th2etaL] = Calc_mu(y, seg, 1, 0, t, 1, 'Local');
    p1 = min(th2etaL(th) ./ qL);  % this is "p1" anchor which is the smallest local bound
    BMax(n) = Calc_B(y, seg, 1, 0, t, p1, th);  % "p1" is enough for comparing targets, since B(p1) is strongly correlated to B(pMax)
end
nBest = argmax(BMax);  % the best blob
[yMaxKb, xMaxKb] = ind2sub(sz(1:2),kBigKb(nBest));  % coodinates of the best pixel

% shows Figure 17(a) - plot the selected pixel on "Kb>2" map
thK = 2;
figure; set(gcf,'Name','high Kbs');
imshow(Kb2d > thK);
fprintf('coverage = %g%%\n', sum(Kb(:) > thK)/numel(Kb)*100)
hold on; plot(xMaxKb, yMaxKb, 'yx', 'MarkerSize', 9, 'LineWidth', 1)

% shows Figure 17(b) - plot on segmentation map (cropped)
cropRoi = [171,227,158,248]; % y1,y2,x1,x2
figure; imagesc(seg(cropRoi(1):cropRoi(2),cropRoi(3):cropRoi(4)))
hold on; plot(xMaxKb-cropRoi(3)+1, yMaxKb-cropRoi(1)+1, 'rx', 'MarkerSize', 9, 'LineWidth', 1)
set(gcf,'position',[1129,592,425,249],'colormap',parula(5)); colorbar;

% shows Figure 17(c) - plot on a PCA visualization (cropped)
figure; imshow(cube2rgb(data(cropRoi(1):cropRoi(2),cropRoi(3):cropRoi(4),:)).^2);
set(gcf,'position',[549,497,557,440]); hold on;
plot(xMaxKb-cropRoi(3)+1, yMaxKb-cropRoi(1)+1, 'yx', 'MarkerSize', 9, 'LineWidth', 1)
