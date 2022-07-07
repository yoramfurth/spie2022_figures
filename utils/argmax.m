function k = argmax(x)
%ARGMAX finds the argument of the largest element in "x".
%
%Description: 
%    This function returns the second output of max(x) built-in function, that is the index/es of the 
%    largest element/s. If "x" is a vector only one output will be calculated, corresponding to the 
%    largest element in X. For matrices, it calculates a row vector containing indexes of the maximal 
%    element from each column. For N-D arrays, it operates similarely along the first non-singleton 
%    dimension.
%
%Inputs: 
% 	 x - an N-D array.
% 
%Outputs: 
% 	 k - indexes of the maximal value/es.
% 
%Examples: 
%    k = ARGMAX([4,1,2,7,3]);  % returns 4
%    k = ARGMAX(magic(5));  % returns [2 1 5 4 3]
% 
%See also MAX, ARGMIN
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

k = argmin(-x);