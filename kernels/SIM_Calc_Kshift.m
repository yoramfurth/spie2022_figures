function Kshift = SIM_Calc_Kshift(K, theta, ar)
%SIM_CALC_KSHIFT calculate angular impact of the inhomogeneity given simulation generating factors.
%
%Description: 
% 	 Kshift is the angular component of "Kb", the directional inhomogeneity impact, given by a local
%    to global ratio of the NSMF expectations. This function calculates it analytically directly from
%    the simulation generating factors. Details are available in the thesis report Section 5.4.2, 
%    regarding Equation 5.38.
% 
%Inputs: 
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0.
%    ar - Aspect-ratio of the ellipsoid major-to-minor axis of a covariance matrix.
% 
%Outputs: 
% 	 Kshift - Directional non-stationarity impact of type shift.
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
ang = theta/2;
dKsp = 0;
Ksplit = K;  % the scaling component of "Kb"
ca2 = (cosd(ang)).^2; 
sa2 = (sind(ang)).^2;

% bias to to Ksplit
if (exist('Ksplit','var') && any(Ksplit~=1))
    Ksplit2 = Ksplit.^2;
    KspNorm = ((Ksplit2 - 1) ./ (Ksplit2 + 1)).^2;  % a normalized version of Ksplit
    dKsp = KspNorm .* (ar.^2 - 2 + ar.^-2) .* ca2 .* sa2;  % an internal shift due to Ksplit
end

% Applies Kshift formula
KG = ca2 + ar.^-2 .* sa2; % main part of XS_G
K0 = ca2 + ar.^2 .* sa2;  % main part of XS_0
Kshift =  (K0 .* KG - dKsp).^0.5;

