function [V,d]=eig_sorted(A)
%EIG_SORTED calculates sorted eigenvalues of a normal matrix, like the covariance matrix. 
%
%Description: 
%    This function produces a column vector "V" containing the eigenvalues of a square matrix "A",
%    and sorted by the eigenvalues "d", in descending order. "A" is assumed to be a normal matrix 
%    (e.g. a covariance matrix), and therefore is ensured to have  a real valid decomposition.
%    Note that this function uses EIG() which is more accurate SVD(). But since it is not guaranteed
%    to get a specific order, we reorganize the provided data by the eigenvalues. 
% 
%Inputs: 
% 	 A - A normal matrix (e.g. a covariance matrix).
% 
%Outputs: 
% 	 V - Eigenvectors, sorted by descending "d".
%    d - A vector of the respective eigenvalues.
% 
%Example: 
%    phi = cov(x);
%    [V,d] = eig_sorted(phi);  % get a transformed data
% 
%See also SVD, EIG
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% extracts the eigenvectors and eigenvalues 
[V,D]=eig(A); 

% organize in descending order
[d,kSort]=sort(diag(D),'descend');
V = V(:,kSort);  
