function k = argmin(x)
%ARGMIN finds the argument of the smallest element in "x".
%
%Description: 
%    This function returns the second output of min(x) built-in function, that is the index/es of the 
%    smallest element/s. If "x" is a vector only one output will be calculated, corresponding to the 
%    smallest element in X. For matrices, it calculates a row vector containing indexes of the minimum 
%    element from each column. For N-D arrays, it operates similarely along the first non-singleton 
%    dimension.
%
%Inputs: 
% 	 x - An N-D array.
% 
%Outputs: 
% 	 k - Indexes of the minimal value/es.
% 
%Examples: 
%    k = ARGMIN([4,1,2,7,3]);  % returns 2
%    k = ARGMIN(magic(5));  % returns [3 2 1 5 4]
% 
%See also MIN
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

[~,k] = min(x);