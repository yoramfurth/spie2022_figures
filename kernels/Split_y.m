function y = Split_y(y, seg, K, theta)
%SPLIT_Y splits data by simulation generating factors.
%
%Description: 
%    For real data simulations we use a fully stationary data as a baseline which have 2 simple segments, 
%    e.g. SIC data-cube. This data is then split to form any inhomogeneity as desired. This function 
%    implements this split given two factors: scaling (K), rotation (theta). Scaling is done by 
%    multiplying the residual noise (y) of one segment by K, and then normalizing the whole data by 
%    sqrt(2/(K^2+1)) to maintain the same global covariance. Rotation is done by rotating each segment by
%    opposite angles of theta/2, along the major-to-minor global eigenvectors' plane. Note that the 
%    implementation assumes that both segments have the same area size. More details about this rotation
%    are available in the thesis report Section 5.3.3. 
% 
%Inputs: 
% 	 y - The residual data to process, fully stationary, with 2 segments.
%    seg - Segmentation map with 2 labels.
%    K - Total scaling of the two segments. Null is 1. 
%    theta - Total rotation angle of the two segments, in degrees. Null is 0. 
% 
%Outputs: 
% 	 y - The residual data after splitting the two segments as to get different covariances. 
% 
%Example: 
%    seg = kmeans_robust(data, 2);  % robustly segments data into 2 clusters
%    y = Calc_y(data, seg, 1, 1);  % calculates the residual noise per-segment, and aligns the stationarity
%    y = Split_y(y, seg, 4, 20);  % splits into 2 distinct segments with x4 scaling, and a 20deg rotation gap 
% 
%See also CALC_Y
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% splits by angle
y = Split_y_byAngle(y, seg, theta, K); % must be first, as implementation assumes that both of segments have the same covariance matrix 

% splits by scaling
y = Split_y_byScaling(y, seg, K);


function y = Split_y_byScaling(y, seg, K)
%SPLIT_Y_BYSCALING multiplies the residual noise (y) of one segment by K, and then normalizes the whole data.

% init.
if (K==1) % time saving
    return;
end

% Multiplies segment #2 by K
k = (seg(:)==2);
y(k,:) = y(k,:)*K; % splits the segment 2 by factor of K

% Normalizes the whole data to maintain the same global covariance (phiG)
y = y * sqrt(2/(K^2+1));  % assumes equal areas in both segments (N1=N2)


function y = Split_y_byAngle(y, seg, theta, K)
%SPLIT_Y_BYANGLE rotates each segment by theta/2 along the major-to-minor plane.
% Assumption: both segments have identical covariance matrix

% init.
if (theta==0) % time saving
    return;
end

% Determines the rotation plane
phiG = cov(y);  % estimates the global covariance matrix from the residual data 
[E,d] = eig_sorted(phiG);  % eigenvectors decomposition 
ax1 = 1;  % index of the major eigenvctor, axis #1 of the rotation plane.
ax2 = numel(d);  % index of the minor eigenvctor, axis #2 of the rotation plane.

% Set opposite angles as to rotate symmetrically and get a gap of "theta" degrees 
theta1 = -theta/2;
theta2 = theta/2;

% Rotates each segment
kSeg = (seg(:)==1);
y(kSeg,:) = rotateData(y(kSeg,:), theta1, E, ax1, ax2);

kSeg = (seg(:)==2);
y(kSeg,:) = rotateData(y(kSeg,:), theta2, E, ax1, ax2);

% Normalizes the whole data to maintain the same global covariance (phiG)
ar2 = d(ax1) / d(ax2);  % AR^2, given by major-to-minor aspect-ratio
dKG = ((K^2 - 1)/(K^2 + 1))^2 * (((1 + 1/ar2)^2) / (cosd(theta2)^-2 + sind(theta2)^-2 / ar2));
Scl = (cosd(theta2)^2 + sind(theta2)^2 / ar2 - dKG)^-0.5;  % see thesis report Equation 5.38
y = y * Scl;


function y = rotateData(y, ang, E, ax1, ax2)
%ROTATEDATA rotates data along a requested plane. 
%
%Description:
%    The equested plane is detemined by 2 eigenvectors (ax1,ax2) from the eigenvectors matrix (E).
%    The rotation is done by determining a rotation matrix for those pair of vectors only, upon
%    an indentity matrix. See details in thesis report Section 5.3.3.
%
%Inputs: 
%    y - Data to rotate.
%    ang - Angle to rotate, in degrees.
%    E - Initial basis (aurthonormal).
%    ax1, ax2 - Indices of axes to rotate.
% 
%Outputs: 
%    y - Rotated data.
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
if (ang==0)
    return;
end

% rotates the requested axes
e1 = E(:,ax1);
e2 = E(:,ax2);  % minor axis
ea1 = cosd(ang) * e1 + sind(ang) * e2;
ea2 = -sind(ang) * e1 + cosd(ang) * e2;

% creates a new basis 
Ea = E;
Ea(:,ax1) = ea1;
Ea(:,ax2) = ea2;

% applies the rotation
y = y * E * Ea';

