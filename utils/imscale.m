function J = imscale(I, inRange, outRange)
%IMSCALE scales image intensity values. 
%
%Description: 
%   J = IMSCALE(I,[LOW_IN; HIGH_IN],[LOW_OUT; HIGH_OUT]) maps the values
%   in data I to new values in J such that values between LOW_IN and HIGH_IN 
%   are mapped to values between LOW_OUT and HIGH_OUT. 
%
%   Note that this function is Similar to imadjust(), but is designed for  
%   general linear transforms, on any arbitrary arrays. 
%
%Inputs: 
% 	 I - Input data, of any dimensionality. 
%    inRange - Inputs range (2 scalars).
%    outRange - Output range (2 scalars).
% 
%Outputs: 
% 	 J - The adjusted data.
% 
%Example: 
%    I = double(imread('pout.tif'))/255;
%    K = imscale(I,[0.3 0.7],[0 1]);  % maps [0.3,0.7] to [0,1]
%    figure, imshow(uint8(K))
%
%See also IMADJUST
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

sf = diff(outRange) ./ diff(inRange); % scale factor
J = (double(I) - inRange(1)) .* sf + outRange(1);

