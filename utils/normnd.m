function z = normnd(x, dim)
%NORMND 2-norm of the a requested dimension.
%
%Description: 
%    This function calculate the 2-norm of a data x. But unlike NORM built-in function, it does 
%    it on any requested dimension. Furthermore, in the default case, it does it on the last 
%    dimension of the data, automatically.
% 
%Inputs: 
% 	 x - Any multi-dimensional data. 
%    dim - The dimension/s on which to perform the norm [Optional, default is last dimension].
% 
%Outputs: 
% 	 z - The resulted norm. The "dim" dimension are collapsed. 
% 
%Example: 
%    B = normnd(magic(3),1);  % [9.43, 10.34, 9.43]
% 
%See also NORM, SUM
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

if (nargin<2)
    dim = ndims2(x);
end
y = x.^2;
y = sum(y,dim);
z = sqrt(y);


function dim = ndims2(x)
%NDIMS2 Gives the last valid dimension.
%
%Description: 
%    This function gives the last valid dimension, but unlike the built-in ndims(),
%    here, 1D vectors are also supported.
%
%Inputs: 
% 	 x - Any multi-dimensional data. 
% 
%Outputs: 
% 	 dim - The last valid dimension.
% 
%Example: 
%    dim = ndims2(magic(3));  % 2
% 
%See also NDIMS
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

dim = max([1, find(size(x)>1,1,'last')]);
