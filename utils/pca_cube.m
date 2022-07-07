function B = pca_cube(A, whiten)
%PCA_CUBE applies PCA transform on a multi-layer image.
%
%Description: 
% 	 This function applies PCA transform using an estimated covariance matrix. It is 
%    then decomposed using SVD, sorted by its eigenvalues. The  resulted eigenvectors 
%    are finally used for transforming the data. If "whiten" is selected, then the data 
%    is also normalizes by the eigenvalues (STD of each layer).
% 
%Inputs: 
% 	 A - A data of at least 2 dimensions. The transform is applied on its last dimension.
% 	 whiten - Do whitening. default is 'false' [optional].
% 
%Outputs: 
% 	 B - Transformed data, with the same dimensions as "A".
% 
%Example: 
%    B = pca_cube(A);  % creates a transformed data
% 
%See also PCA, SVD, EIG
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

if (nargin<2)
    whiten = false;
end
sz = size(A);
A = reshape(A,[prod(sz(1:end-1)),sz(end)]); % ensure 2D
covMat = cov(A);
[V,lam] = eig_sorted(covMat);
if (~whiten)
    B = A * V;
else
    B = A * V * diag(1./sqrt(lam));
end
B = reshape(B,sz);
