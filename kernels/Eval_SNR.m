function snr = Eval_SNR(data, y)
%EVAL_SNR evaluates the average SNR of a data-cube.
%
%Description: 
%    This function evaluates the average SNR directly from the data-cube.
% 	 The average SNR refers here to ratio between an average data sample, and 
%    a typical noise sample. These two values are estimated from the data. The 
% 	 former is by averaging the magnitude of all the pixels in the data-cube. 
%    The latter is similar but on the residual noise of that cube, that is on
%    the "y=x-m" cube (see readme.txt). 
% 
%Inputs: 
% 	 data - A data-cube. 
%    y - The residual noise of each pixel.
% 
%Outputs: 
% 	 snr - The average SNR of the data-cube.
% 
%Example: 
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    avgSNR = Eval_SNR(data, y);  % evaluates the global SNR
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

averageDataSample = mean2(normnd(data));  % mean magnitude of all the pixels
averageNoiseSample = mean(normnd(y));  % mean magnitude of all the residual samples
snr = averageDataSample / averageNoiseSample;  
