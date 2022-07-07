function mm = minmax2(x)
%MINMAX2 claulates [min,max] of the given data.
%
%Description: 
%   This function is a combination of min(),max(). It finds the bounds of
%   an array regardless of its type or dimensionality (unlike minmax). 
%
%Inputs: 
% 	 x - a data of any kind.
% 
%Outputs: 
%    mm - the minimal and maximal value.
%    
%Example: 
%        X = [2 8 4 ;7 3 9];
%        mm = minmax2(X)  % gives [2 9]
%    
%See also MIN, MAX, MINMAX
%    

%    Copyright 2012-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

mm = [min(x(:)), max(x(:))];
