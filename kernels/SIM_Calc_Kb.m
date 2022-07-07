function Kb = SIM_Calc_Kb(phiG, phiS, t)
%SIM_CALC_KB calculates the directional inhomogeneity impact.
%
%Description: 
% 	 Kb is a new metric that estimates the impact of inhomogeneity along the disctiminant axis.
%    This measure is based on comparing the maximal local expectation to the global expectation. 
%    Local expectation refers to with-target NSMF expectation per segment. Global expectation
%    refers to with-target NGMF expectation.
% 
%Inputs: 
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
% 	 t - A matrix representing an array of vectors representing targets spectrum. 
% 
%Outputs: 
% 	 Kb - Directional inhomogeneity impact.
% 
%Example: 
%    phiS = {[1 0; 0 0.25]; [0.25 0; 0 1]}
%    phiG = 0.5 * phiS{1} + 0.5 * phiS{2};
%    t = [20 90]';
%    Kb = SIM_Calc_Kb(phiG, phiS, t);
%
%See also SIM_CALC_MU
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

qG = Eval_q(t, phiG);  % evaluates the full-target global expectation  
qL = cellfun(@(phi)Eval_q(t,phi), phiS, 'UniformOutput', false);  % evaluates the full-target local expectation per each segment
q1 = max(vertcat(qL{:}));  % selects the maximal local expectation
Kb = q1 ./ qG;  % the direction inhomogeneity measure

