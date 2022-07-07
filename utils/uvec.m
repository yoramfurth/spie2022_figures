function uvec = uvec(vec)
%UVEC calculates the unit vector of a vector by normalizing it by its magnitude.
%
%Description: 
%    A unit vector is a vector of magnitude 1. This can be achieved by normalizing 
%    a vector by its magnitude. This function does this at once on a set of vectors.
%    The function assumes that the vectors lies on the last dimension, the vetors 
%    might be orgnized in any N-D array from the moment that the last dimension is 
%    for the vectors. For example in an M-by-N-by-D data, only the "D" dimension
%    will be normalized. 
% 
%Inputs: 
% 	 vec - A vector, or a set of vectors (vectors are the last dimension).
% 
%Outputs: 
% 	 uvec - Normalized vectors (dimensions are preserved).
% 
%Example: 
%    u = uvec(magic(5));  % unit vectors
%    nrm = normnd(u)';  % all magnitudea gets 1.0
% 
%See also NORM, NORMND
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

uvec = vec ./ normnd(vec);
