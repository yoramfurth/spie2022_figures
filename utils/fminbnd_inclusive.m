function [xf, fval] = fminbnd_inclusive(fun, ax, bx, varargin)
%FMINBND_INCLUSIVE Single-variable bounded nonlinear function minimization, bounds inclusive. 
%
%Description: 
%    This function minimizes fun() in [ax,bx] bounds, but unlike the builtin fminbnd(), it also
%    checks the bounds themselves.  
% 
%Inputs: 
% 	 fun - A function to minimize, which receives one scalar argument.
%    ax - The lower bound of the scanned interval. 
%    bx - The upper bound of the scanned interval. 
% 
%Outputs: 
% 	 xf - The argument that give the minimal value, bounds inclusive. 
%    fval - The minimal value, bounds inclusive. 
% 
%See also FMINBND
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

varsOut = cell(1,max(1,nargout));
[xf, fval] = fminbnd(fun, ax, bx, varargin{:});

xf_inc = [ax, xf, bx];
[fval_inc,k] = min(arrayfun(fun, xf_inc));
if (k~=2)  
	xf = xf_inc(k);
    fval = fval_inc;
% else, remains with the result of fminbnd()
end
