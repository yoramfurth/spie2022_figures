function [phi1, phi2] = SIM_Split_phi(phiG, K, theta)
%SIM_SPLIT_PHI splits the global covariance matrix by simulation generating factors.
%
%Description: 
%    For synthetic simulations we use spectral factors that represents the second moment the a data 
%    residual noise. For simulating inhomogeneity we starts with a basic global covariance matrix (phiG), 
%    and "split" it in an equivalent way to real data in SPLIT_Y. In the scaling case (K) this is done by 
%    determining on of the segments by phiG*K^2, and then normalizing both segments by 2/(K^2+1) to 
%    maintain the same global covariance. Rotation is done by rotating each segment by opposite angles of 
%    theta/2. Note the square relationships with the equivalent operations in SPLIT_Y, which occurs  
%    because here we work on the covariances which have a square relationship with the data. Note also
%    that the implementation assumes that both segments have the same area size.
% 
%Inputs: 
% 	 phiG - The global covariance matrix.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 
%Outputs: 
% 	 phi1, phi2 - The covariance matrix of segment #1 and #2, respectively.
% 
%Example: 
%    AR = 120; % major-to-minor aspect-ratio of the covariance matrix
%    phiG = SIM_Get_phi(AR);
%    [phi1, phi2] = SIM_Split_phi(phiG, 4, 20);  % determines 2 distinct segments with x4 scaling, and a 20deg rotation gap 
% 
%See also SPLIT_Y
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
if (~iscell(phiG))
    [phi1, phi2] = deal(phiG);
else
    [phi1, phi2] = deal(phiG{:});
end

% splits by angle
[phi1, phi2] = SIM_Split_phi_byAngle(phi1, phi2, theta, K);

% splits by scaling
[phi1, phi2] = SIM_Split_phi_byScaling(phi1, phi2, K);


function [phi1, phi2] = SIM_Split_phi_byScaling(phi1, phi2, K)
%SIM_SPLIT_PHI_BYSCALING multiplies "phi" of one segment by K, and then normalizes the whole data.

% init.
if (K==1)
    return;
end

% Multiplies segment #2 by K
phi2 = phi2 * K^2;

% Normalizes the whole data to maintain the same global covariance (phiG)
Scl = 2/(K^2+1);  % assumes equal areas in both segments (N1=N2)
phi1 = phi1 * Scl;
phi2 = phi2 * Scl;


function [phi1, phi2] = SIM_Split_phi_byAngle(phi1, phi2, theta, K)
%SPLIT_Y_BYANGLE rotates each segment by theta/2.

% init.
if (theta==0)
    return;
end

% Set opposite angles as to rotate symmetrically and get a gap of "theta" degrees 
theta1 = -theta/2;
theta2 = theta/2;

% Set normalization factors to the whole data as to maintain the same global covariance (phiG)
dKG = @(ar2) ((K^2 - 1)/(K^2 + 1))^2 * (((1 + 1/ar2)^2) / (cosd(theta2)^-2 + sind(theta2)^-2 / ar2));
fScl = @(ar2) 1/(cosd(theta2)^2 + sind(theta2)^2 / ar2 - dKG(ar2));  % see thesis report Equation 5.38

[~,R,~] = svd(phi1);  % svd gets descending order
ar2 = R(1,1)/R(end,end);  % AR^2, given by major-to-minor aspect-ratio
Scl(1) = fScl(ar2);

[~,R,~] = svd(phi2); % svd gets descending order
ar2 = R(1,1)/R(end,end);  % AR^2, given by major-to-minor aspect-ratio
Scl(2) = fScl(ar2);

% Applies scaling and rotation on the covariance matrixes
phi1 = RotateCov2D(phi1 * Scl(1), theta1);
phi2 = RotateCov2D(phi2 * Scl(2), theta2);


function C2 = RotateCov2D(C1,alpha)
%ROTATECOV2D rotates a 2D covariance matrix.
%
%Description: 
%    The function builds a rotation matrix and multiplies the input matrix.
%    The PCA components are proportional to S=C^0.5. Therefore the rotation
%    matrix has to be applied twice, one for each S. Note that the rotation 
%    is CCW for xy order, or CW for ij order.
%
%Inputs: 
%    C1 - input covariance matrix
%    alpha - angle to rotate by [deg]
% 
%Outputs: 
%    C2 - output covariance matrix
%   

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

R = [cosd(alpha), -sind(alpha);...
    sind(alpha), cosd(alpha)];
C2 = R * C1 * R';
