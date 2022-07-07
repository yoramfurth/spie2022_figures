function [A, Anom, dAnom] = SIM_Eval_A(cN, cT, fN, eta)
%SIM_EVAL_A evaluates how much relative superior detection is achieved over a range of thresholds.
%
%Description: 
%    This function evaluates a score for the ROC-curve within a range of Pfa values, as introduces in EVAL_A.  
%    This function calculates this analytically, and returns handles to functions of the target's power (p). 
%    More details are available in the paper Section 1.4, and in the thesis report Sections 2.6 and 6.3.
% 
%Inputs: 
% 	 cN(z) - CDF of the scores of NSMF with no target.
% 	 cT(z,p) - CDF of the scores of NSMF with target.
% 	 fN(z) - PDF of the scores of NSMF with no target.
% 	 eta - Decision threshold on the scores axis. In case of array, any value is processed. 
% 
%Outputs: 
% 	 A(p) - Evaluation of the detection algorithm, A(th), for any of the requested "th".
% 	 Anom(p) - The numerator of A(th). Usable for calculating B(th) accurately. 
% 	 dAnom(p) - Optional. Derivative of Anom. Usable for L'hopital's approximation.
% 
%Example: 
%    p = 0.02;  % target's power to test
%    th = 0.01;  % desired number of false alarms
%    phi1 = [1 0; 0 0.25];
%    phi2 = [0.25 0; 0 1];
%    t = [20 90]';
%    [fN,~,~,cN,cT,icN] = SIM_Calc_PDF(phi1, phi2, t, 0.5, 'Local'); 
%    eta = -icN(th);
%    A = SIM_Eval_A(cN, cT, fN, eta)
%    fprintf('A(th)=%g for a threshold of th=%g, and a target's power of p=%g\n', A(p), th, p);
% 
%See also SIM_CALC_PDF, EVAL_A
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% repeats for any "p" and any "eta" 
kAll = @(p) 1 : max(numel(eta), numel(p));  % indexes of all elements
A = @(p) arrayfun(@(k) Build_A(cN, cT, fN, eta(min(k,end)), p(min(k,end))), kAll(p));
Anom = @(p) arrayfun(@(k) Build_Anom(cN, cT, fN, eta(min(k,end)), p(min(k,end))), kAll(p));
dAnom = @(p) arrayfun(@(k) Build_dAnom(cN, cT, eta(min(k,end)), p(min(k,end))), kAll(p));


function A = Build_A(cN, cT, fN, eta, p)
%BUILD_A evaluates "A" for a given thrshold (eta) and a given target's power (p).
funDenom = @(z)cN(z).*fN(z);
Adenom = integral(funDenom, eta, inf);

% check if is stable numerically
if (Adenom>0.1) % high denominator implies stable - evaluates "A" directly
    Anom = Build_Anom(cN, cT, fN, eta, p);
    A = Anom ./ Adenom;
    
else % numeric solution for low denominators - evaluates "A" through the complementary value "1-A"
    funNom = @(z)cT(z,p).*fN(z);
    nom = integral(funNom, eta, inf);
    Acomp = nom ./ Adenom;  % a complementatry value
    dAcomp = cT(eta,p)/cN(eta);  % derivative, for L'Hôpital's approximation at 0/0 limit
    if (isnan(Acomp) || dAcomp<1e-3)
        Acomp = dAcomp;
    end
    A = 1 - Acomp;
end

assert(A>=0 && A<=1 && isfinite(A), 'invalid A');


function Anom = Build_Anom(cN, cT, fN, eta, p)
%BUILD_ANOM evaluates the numerator of "A".
funNom = @(z)(cN(z)-cT(z,p)).*fN(z);
Anom = integral(funNom, eta, inf);


function dAnom = Build_dAnom(cN, cT, eta, p)
%BUILD_ANOM evaluates the derivative of the numerator of "A".
dAnom = cN(eta)-cT(eta, p);  % Derivative of Anom, for L'Hôpital's approximation. Assuming Pfa(end-1)==0, and Pd(end-1)==Pd(end),  


