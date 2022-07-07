function segKMeans = kmeans_robust(X, k)
%KMEANS_ROBUST robust K-means clustering of a data-cube. 
%
%Description: 
%    This function clusters a given "X" data-cube based on K-means algorithm. It 
%    partitions the points in the H-by-W-by-D data-cube "X" into "k" clusters. This
%    partition minimizes the sum, over all clusters, of the within-cluster sums of 
%    point-to-cluster-centroid distances. 
% 
%    This function combines several features to provide this robustly and efficiently.
%    The first key point for this is initializing the clusters with the robust Otsu 
%    thresholding. Unlike the built-in GRAYTHRESH this version supports multi-clusters, 
%    in a way that was proposed in the original Otsu's paper. The second key point is a 
%    simplified implementation of K-means that is faster by many orders, comparing to 
%    the built-in KMEANS function. 
% 
%Inputs: 
% 	 X - A 3-D data-cube, of size [Height, Width, Depth].
% 	 k - The number of requested clusters.
% 
%Outputs: 
% 	 segKMeans - Transformed image.
% 
%Example: 
%    seg = kmeans_robust(data, 5);  % robustly segments data into 5 clusters 
% 
%See also GRAYTHRESH, KMEANS
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).
%
%    graythreshN() - Copyright 2007 Damien Garcia (see details below) 
%    kmeans_fast() - Copyright 2017 Mo Chen (see details below) 

% Initial clustering using Otsu multi-thresholding on mean image 
imMean = mean(X,3);
segOtsu = graythreshN(imMean,k);

% clustering using kmeans 
sz = size(X);
X = reshape(X,[prod(sz(1:2)),sz(3)]); % column-stack
segKMeans = kmeans_fast(X.', segOtsu(:)');
segKMeans = reshape(segKMeans,sz(1:2));
end


function [IDX,sep] = graythreshN(I,n)
%GRAYTHRESHN Global image clustering using Otsu's method.
%
%Description: 
%    This function clusters the image I into "n" classes by means of Otsu's
%    N-thresholding method. GRAYTHRESHN returns an array IDX containing the 
%    cluster indices of each point. Non-finite pixels are assigned to zeros.
%
%    Note that the thresholds generally become less credible as the number of
%    classes (N) to be separated increases (see Otsu's paper for more details).
%
%Inputs: 
%    I - The input image (2D).
%    n - Number of requested clusters. 
%
%Outputs: 
%    IDX - An array containing the cluster indices (from 1 to n) of each point.
%    sep - The value of the separability criterion within the range [0 1].
%
%Example: 
%    I = double(imread('pout.tif'))/255;
%    IDX = graythreshN(I,5);
%    figure; imagesc(IDX);
%    colormap(gray)
%
%Reference: 
%    Otsu N, <a href="matlab:web('http://dx.doi.org/doi:10.1109/TSMC.1979.4310076')">A Threshold Selection Method from Gray-Level Histograms</a>,
%    IEEE Trans. Syst. Man Cybern. 9:62-66;1979
%
%See also GRAYTHRESH, IM2BW
%

%    Copyright 2007/08 Damien Garcia, as "otsu.m" (see license in license_Damien_Garcia.txt)
%    https://www.mathworks.com/matlabcentral/fileexchange/26532-image-segmentation-using-otsu-thresholding
%    
%    Updated by Damien Garcia, 2010/03
%    Updated by Yoram Furth, 2017, and renamed to "graythreshN"

%% Checking n (number of classes)
assert(ismatrix(I), 'The input must be a 2-D array.');
n = round(n);
assert(n>1, 'should have at least 2 cluster');
if n > 255
    n = 255;
    warning('MATLAB:otsu:TooHighN','n is too high. n value has been changed to 255.')
end

%% Convert to 256 levels
I = single(I);
I = I-min(I(:));
I = round(I/max(I(:))*255);

%% Probability distribution
unI = sort(unique(I));
nbins = min(length(unI),256);
if nbins==n
    IDX = ones(size(I));
    for i = 1:n, IDX(I==unI(i)) = i; end
    sep = 1;
    return
elseif nbins<n
    IDX = NaN(size(I));
    sep = 0;
    return
elseif nbins<256
    [histo,pixval] = hist(I(:),unI);
else
    [histo,pixval] = hist(I(:),256);
end
P = histo/sum(histo);
clear unI

%% Zeroth- and first-order cumulative moments
w = cumsum(P);
mu = cumsum((1:nbins).*P);

%% Maximal sigmaB^2 and clustered image
if n==2
    sigma2B =...
        (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
    [maxsig,k] = max(sigma2B);
    
    % clustered image
    IDX = ones(size(I));
    IDX(I>pixval(k+1)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
elseif n==3
    w0 = w;
    w2 = fliplr(cumsum(fliplr(P)));
    [w0,w2] = ndgrid(w0,w2);
    
    mu0 = mu./w;
    mu2 = fliplr(cumsum(fliplr((1:nbins).*P))./cumsum(fliplr(P)));
    [mu0,mu2] = ndgrid(mu0,mu2);
    
    w1 = 1-w0-w2;
    w1(w1<=0) = NaN;
    
    sigma2B =...
        w0.*(mu0-mu(end)).^2 + w2.*(mu2-mu(end)).^2 +...
        (w0.*(mu0-mu(end)) + w2.*(mu2-mu(end))).^2./w1;
    sigma2B(isnan(sigma2B)) = 0; % zeroing if k1 >= k2
    
    [maxsig,k] = max(sigma2B(:));
    [k1,k2] = ind2sub([nbins nbins],k);
    
    % clustered image
    IDX = ones(size(I))*3;
    IDX(I<=pixval(k1)) = 1;
    IDX(I>pixval(k1) & I<=pixval(k2)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
else
    k0 = linspace(0,1,n+1); k0 = k0(2:n);
    [k,y] = fminsearch(@sig_func,k0,optimset('TolX',1));
    k = round(k*(nbins-1)+1);
    
    % clustered image
    IDX = ones(size(I))*n;
    IDX(I<=pixval(k(1))) = 1;
    for i = 1:n-2
        IDX(I>pixval(k(i)) & I<=pixval(k(i+1))) = i+1;
    end
    
    % separability criterion
    sep = 1-y;
    
end

IDX(~isfinite(I)) = 0;

%% Function to be minimized if n>=4
    function y = sig_func(k)
        
        muT = sum((1:nbins).*P);
        sigma2T = sum(((1:nbins)-muT).^2.*P);
        
        k = round(k*(nbins-1)+1);
        k = sort(k);
        if any(k<1 | k>nbins), y = 1; return, end
        
        k = [0 k nbins];
        sigma2B = 0;
        for j = 1:n
            wj = sum(P(k(j)+1:k(j+1)));
            if wj==0, y = 1; return, end
            muj = sum((k(j)+1:k(j+1)).*P(k(j)+1:k(j+1)))/wj;
            sigma2B = sigma2B + wj*(muj-muT)^2;
        end
        y = 1-sigma2B/sigma2T; % within the range [0 1]
        
    end

end


function [label, dstEnergy, dstMeans] = kmeans_fast(X, k_or_labels)
%KMEANS_FAST Perform simplified fast k-means clustering.
%
%Description: 
%    This function partitions the points in the N-by-P data matrix X
%    into K clusters. This partition minimizes the sum, over all clusters, of
%    the within-cluster sums of point-to-cluster-centroid distances. Rows of X
%    correspond to points, columns correspond to variables.
%
%Inputs: 
%    X - d-by-n data matrix
%    k_or_labels - option 1: number of clusters ("K"), option 2: inital labels (1-by-n vector)
% 
%Outputs: 
%    label - 1-by-n cluster label
%    dstEnergy - optimization target's value
%    dstMeans - trained model structure
%
%See also KMEANS
% 

%    Copyright 2017 Mo Chen (sth4nth@gmail.com), as "kmeans.m" (see license in license_Mo_Chen.txt)
%    https://www.mathworks.com/matlabcentral/fileexchange/24616-kmeans-clustering
%    
%    Updated by Yoram Furth, 2017, and renamed to "kmeans_fast"

n = size(X,2);
if numel(k_or_labels)==1
    k = k_or_labels;
    label = ceil(k*rand(1,n));
elseif numel(k_or_labels)==n
    label = k_or_labels;
end

last = 0;
while any(label ~= last)
    [u,~,label(:)] = unique(label);  % removes empty clusters
    k = numel(u);
    E = sparse(1:n,label,1,n,k,n);  % transforms label into indicator matrix
    m = X*(E*spdiags(1./sum(E,1)',0,k,k));  % computes centers 
    last = label;
    [val,label] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);  % assigns labels
end
dstEnergy = dot(X(:),X(:))-2*sum(val); 
dstMeans = m;
end
