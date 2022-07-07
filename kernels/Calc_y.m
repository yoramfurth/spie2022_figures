function y = Calc_y(data, seg, segmentedMean, forceStationary)
%CALC_Y calculates the residual noise of each pixel in the dataset. 
%
%Description: 
% 	 Calculating the residual noise is done basically by subtracting each pixel local average from itself, 
%    as introduces in SUBTRACT_LOCAL_AVERAGE. However one might want to align artificially the local  
%    covariance to each others, as to form a stationary cube (assuming that the data within each segments 
%    is already stationary). This is the case, for example, for SIC data-cube which is used a baseline
%    for simulations, assuming that it is fully stationary (in the second moment sense). 
%
%    Forcing stationarity is done here using whitening techniques. For each segment is first whitened 
%    by its local covariance matrix (phiS), and the unwhiten by the global covariance matrix (phiG). 
%    More details appear in the thesis report, Section 5.3.3.
% 
%Inputs: 
% 	 data - A data-cube. 
%    seg - Segmentation map.
%    segmentedMean - If yes, considers the segmentation while subtracting the local mean. 
%    forceStationary - If yes, align the data of the segment as to have identical covariance matrix. 
% 
%Outputs: 
% 	 y - The residual data to process, basically a subtraction of the local average from each pixel (x-m).
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    y = Calc_y(data, seg, 1, 1);  % calculates the residual noise per-segment, and aligns the stationarity.
% 
%See also SUBTRACT_LOCAL_AVERAGE
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).


% init.
if (ismatrix(data)) % by convention, a 2D data is for "y", whereas original "data" comes in 3D matrix 
    y = data; 
    return;  % "y" has already been calculated
end

% subtract local average
if (~segmentedMean)
    y = Subtract_Local_Average(data);
else
    y = Subtract_Local_Average(data, seg);
end

% reshape images => columns
sz = size(data);
y = reshape(y,[sz(1)*sz(2),sz(3)]);

% fix non-stationarity of segments covariances
if (forceStationary)
    y = allignStationarity(y, seg);
end


function y = allignStationarity(y, seg)
% ALLIGNSTATIONARITY forces data to have same covariance on each segment.
%
% The naming convention here is suffix "G" for global measures, and suffix "S" for 
% per-segment measures.
%

% decomposition
phiG = cov(y);
[EG,dG] = eig_sorted(phiG);

% decompose each segment and transform it onto phiG
seg = seg(:);
for s = unique(seg)' 
    k = (seg==s);  % mask of the current segment
    yS = y(k,:);  % residual in this segment
    phiS = cov(yS);  % covariance matrix of this segment
    [ES,dS] = eig_sorted(phiS);  % eigenvectors decomposition
    y(k,:) = yS * ES * diag(sqrt(dG./dS)) * EG';
end
