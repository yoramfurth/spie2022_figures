function [B, AG, AL] = Calc_B(data, seg, K, theta, t, p, th)
%CALC_B calculates B as function of various factors.  
%
%Description:
%    The metric "B" evaluates the segmentation benefit, achieved over a range of thresholds.  
%    CALC_B receives basic spectral factors and calculates "B" from zero, by computing all 
%    modules  up to the required point. This function can test set of values for any relevant 
%    factor. For example, if "p" is an array, then it results "B(p)" for any "p" in the array. 
% 
%Inputs: 
% 	 data (or y) - The data-cube to process (of the residual noise).
%    seg - Segmentation map.
%    K - Total scaling of the two segments. Null is 1. Supports array.
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. Supports array.
% 	 t - A vector representing the spectrum of the target of interest. Supports array.
%    p - The target's power. Supports array.
%    th - The decision thresholds in terms of Pfa. Supports array.
% 
%Outputs: 
% 	 B - The segmentation benefit. 
% 	 AG, AL - Detection algorithm evaluator in global and local domain, respectively.
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    th = 1./10.^1:0.1:3;  % th from 0.001 to 0.1
%    B = Calc_B(data, seg, 4, 20, t, 1, th);  % test spiting into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
%    figure; loglog(th, B);  % plots B(th) from 0.001 to 0.1
% 
%See also EVAL_B, CALC_A
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

[AG, AGnom, dAGnom] = Calc_A(data, seg, K, theta, t, p, th, 'Global');
[AL, ALnom, dALnom] = Calc_A(data, seg, K, theta, t, p, th, 'Local');
B = Eval_B(AG, AL, AGnom, ALnom, dAGnom, dALnom);
