function q = Eval_q(t, phi)
%EVAL_Q evaluates 
% evaluates mu(p=1), which is the SNR in NSMF.
%
%Description: 
%    Assuming a zero expectation without target, i.e. E(x-m)==0, appending a target (p*t)
%    biases the NSMF and the expectation becomes mu=E(NSMF|T)=p*sqrt(t'*phi^-1*t). The 
%    special case where p=1 is denoted "q", and it refers to the SNR for a full target
%    addition (since NSMF has a unit STD). This measure is a common building block while
%    analyzing NSMF. This function supports also a vactor of many targets. 
% 
%Inputs: 
% 	 t - A matrix representing an array of vectors representing targets' spectrum. 
%    phi - A covariance matrix.
% 
%Outputs: 
% 	 q - NSMF expectation for a full target addition.
% 
%Example: 
%    y = Calc_y(data, [], 0, 0);  % calculates the residual noise (x-m)
%    phiG = cov(y);  % calculates the global covariance  
%    snr = Eval_q(t, phiG);  % evaluates the global SNR
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

q = sqrt(dot(t, phi \ t));  % "dot" supports vector of many targets, i.e. t==[t1,t2,t3,...]

