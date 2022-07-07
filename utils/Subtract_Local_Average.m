function bgEst = Subtract_Local_Average(im, seg)
%SUBTRACT_LOCAL_AVERAGE subtracts local average from each pixel.
%
%Description: 
%    The goal of this function is to subtract the background of a pixel
%    as to bring out sub-pixel targets, and detect them. Here we estimate
%    the background by averaging the 8 neighbors of the pixel.
%    
%    If no segmentation, it's just a simple mean filter, but on the presence
%    of segmentation it is more complex. The problem arises on non-convex borders
%    of the segments. In the one hand we don't want pixel from other segments
%    to participate in the estimation. But in the other hand just stopping 
%    on the borders causes an inherent bias to the global mean. After all, 
%    after subtracting the mean we'd expect to get zero-mean in each segment.
%    
%    This function solves this problem. It is done by a special replicate 
%    which takes into consideration non-convex borders. Basically it's done
%    by a weighted average, but of a new type. Each pixel near the border 
%    gets the average of the nearby pixels within the segment, and also gets 
%    a weight related to their amount. Then both are used for doing the 
%    weighted average on the borders. On the image border it remains exactly
%    like the normal "symmetric" mode. 
%
%Inputs: 
%    im - A data-cube image.
%    seg - Segmentation map of the image [Optional].
%
%Outputs:
%    bgEst - The background estimation for each pixel neighborhood.
%
%Examples:
%    bgWithoutSeg = SUBTRACT_LOCAL_AVERAGE(cubeData);  
%    bgWithSeg = SUBTRACT_LOCAL_AVERAGE(cubeData, mapOfLabels);
%
%See also IMFILTER 
%

%    Copyright 2019-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

if (nargin<2)
    seg = [];
end

% the core algorithm
h = ones(3)/8; h(2,2)=0;
meanEstFun = @(im) imfilter(im,h,'symmetric');  % local mean estimation

% handle segmentation case
[sMax, sAll] = GetNumOfSegments(seg);
if (sMax<=1)
    m = meanEstFun(im);
else
    m = estimateLocalAverage(im, seg, meanEstFun, sAll);
end

% subtracts the estimated local average
bgEst = im - m;


function m = estimateLocalAverage(im, seg, meanEstFun, sAll)
% Calculates the local mean per segment.
% It's done by replicating the perimeter of the segment
% This is like what "symmetric" does, but adapted to non-convex borders. 
% 
Nbands = size(im,3);
sMax = numel(sAll);
m = zeros(size(im),'like',im);
for ks = 1:sMax
    kSeg = (seg==sAll(ks));
    im_seg = replicateSegPerimeter(im, kSeg, meanEstFun);
    mInSeg = meanEstFun(im_seg);
    k = find(repmat(kSeg, [1 1 Nbands]));
    m(k) = mInSeg(k);
end


function im = replicateSegPerimeter(im, kSeg, meanEstFun)
% replicates the perimeter of the segment
im = im.*kSeg;
h_norm = max(eps, meanEstFun(double(kSeg)));  
im = im + meanEstFun(im) ./ h_norm .* (~kSeg);  %  im += (im**h).*(~seg)


function [sMax,sAll]=GetNumOfSegments(seg)
%GETNUMOFSEGMENTS calculate the number of segments (sMax) and their labels (sAll) 
if isempty(seg)
    sAll = [];
    sMax = 1;
else
    sAll = unique(seg(:));
    sMax = numel(sAll);
end
