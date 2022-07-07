function Show_Figure14c(SIC_data, SIC_seg)
%SHOW_FIGURE14C shows Figure 14(c) - Influence of the target's direction in the two-segments case on a representative layout.
%
%Inputs: 
% 	 SIC_data - SIC data-cube. 
% 	 SIC_seg - SIC data-cube's segmentation (2 segments). 
%
%See also CALC_PMAX
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
theta = 90;
SNR = 100; 
th = 0.01;
AR = 8;
seg = SIC_seg;

% sampling ranges
beta_r1 = 1./[1,2,6];  % scaling factor of r1 (the first ellipse's axis)
aAll = -10:0.5:(theta+10);  % target's rotation angle

% extracts minor and major axes
y = Calc_y(SIC_data, seg, 1, 1);
phiG = cov(y);  % global covariance
E = eig_sorted(phiG);  % eigenvectors decomposition
yax1 = E(:,1);  % major axis, defines axis1 of the rotation plane
yax2 = E(:,end);  % minor axis, defines axis2 of the rotation plane

% init. loop
Ntests = max(numel(beta_r1));
Bmax = zeros(Ntests,numel(aAll));

% for each scaling factor calculates Bmax
for n=1:Ntests
    fprintf('test %d/%d', n, Ntests);
	yn = Split_y_as_ellipsoids(y, seg, phiG, AR, [0,theta], {[1 1], [1, beta_r1(n)]});  % extract the data background
    
    for m=1:numel(aAll)
        if (mod(m,20)==0)
            fprintf('.');
        end
		
		% defines the target that corresponds to this angle
		a = aAll(m);
		t = cosd(a)*yax1 + sind(a)*yax2;
        [MF_N, MF_T] = Calc_MF(yn, t, 'Global', seg, 1);
        [qG, ~] = Eval_mu(MF_T, seg, 'Global');
        t = t .* (SNR / qG);  % fix SNR

		% finds the local p-anchors
		pG = th2eta(MF_N(:),th) / SNR;
        [MF_N, MF_T] = Calc_MF(yn, t, 'Local', seg, 1);
        [q0, qL] = Eval_mu(MF_T, seg, 'Local');
        pL = max(eps, th2eta(MF_N(:),th)) ./ qL;
        
        % scans B(p) in the range of p-anchors (similar to Calc_pmax())
        pAnc = [pL(:); pG];  % all p-anchors
		minRng = min(pAnc)/10;
		maxRng = max(pAnc)*10;
		pRng = 10.^(log10(minRng): 0.01 :log10(maxRng));
		B = Calc_B(yn, seg, 1, 0, t, pRng, th);

		% estimates pMax,Bmax (similar to Calc_pmax())
		B0 = imfilter(B, fspecial('gaussian',19,3), 'symmetric');
		Bmax(n,m) = max(B0);  % estimates Bmax. Very similar 
    end
    fprintf('\n');
end
fprintf('done!\n\n');

% init. plots
axMrgP = @(ax)([ax(1:2),10.^(log10(ax(3:4))+diff(log10(ax(3:4)))*0.02*[-1,1])]);
vec2legendcell = @(v, vname) split(strtrim(sprintf([vname,'=1/%.0g '], 1./v)),' ');

% plots Bmax(a)
figure;
hPlot = semilogy(aAll, Bmax);
xlabel('\alpha [deg]'); ylabel('Bmax'); grid on;  % a is denoted "alpha" in the documents
ax = [aAll([1,end]), min(Bmax(:)/1.3), max(Bmax(:))]; axis(axMrgP(ax));
hold on; plot([90 90], ax(3:4), 'k--'); hold on; plot([0 0], ax(3:4), 'k--'); uistack(hPlot,'top');
legend(hPlot, vec2legendcell(beta_r1, '\beta_{a}'),'location','SouthWest');  % r1 is called "a" in the documents


function y = Split_y_as_ellipsoids(y, seg, phiG, AR, rotAngs, scl)
%SPLIT_Y_AS_ELLIPSOIDS splits data as to get two covariances' ellipsoids with the desirable components.
%
%Description: 
%    Using a respective linear transform, its possible to deform any data to get any desirable covariance. 
%    Specifically, this function transforms each segment along the major-to-minor global eigenvectors' plane,
%    using 2D scaling and rotation operations. Along this plane, each covariance gets a desirable aspect-ratio, 
%    a desirable rotation angle, and a desirable extra scaling per axis. Note that the implementation assumes
%    that "y" is fully stationary, i.e. both segments have the same covariance. Note also that the it assumes
%    that both segments have the same area size. More details are available in the thesis report Section 8.2.3,
%    regarding Equation 8.24.
% 
%Inputs: 
% 	 y - The residual data to process, fully stationary, with 2 segments.
%    seg - Segmentation map with 2 labels.
%    phiG - The global covariance matrix. This is the baseline of the transforms. 
%    AR - Desirable aspect-ratios of both covariances (phi1, phi2). 
%    rotAngs - Rotation angle of the each segment, in degrees. Null is 0. 
%    scl - Extra scaling per segment per 2D ellipse's axis. Null is 1. 
%
%Outputs: 
% 	 y - The residual data after splitting the two segments. 
% 
%See also SIM_SPLIT_PHI, SPLIT_Y
%

% init.
[E,d] = eig_sorted(phiG);
L = numel(d); % L, after number of Layers
ax1 = 1;  % index of major axis, which defines axis1 of the rotation plane
ax2 = L; % index of minor axis, which defines axis2 of the rotation plane

% scaling
scalingFactor = AR / sqrt(d(ax1) / d(ax2));
K = sqrt(scalingFactor) .^ [1, -1];  % splits scaling to each axis

% assuming 2 segments
S_all = unique(seg(:));
kSeg1 = (seg(:)==S_all(1));
kSeg2 = (seg(:)==S_all(2));

% transform each segment
y(kSeg1,:) = transformData(y(kSeg1,:), K, rotAngs(1), E, ax1, ax2);
y(kSeg2,:) = transformData(y(kSeg2,:), K, rotAngs(2), E, ax1, ax2);

% post-scaling
y(kSeg1,:) = scalingData(y(kSeg1,:), E, ax1, ax2, scl{1});
y(kSeg2,:) = scalingData(y(kSeg2,:), E, ax1, ax2, scl{2});


function y = transformData(y, K, rotAng, E, ax1, ax2)
%TRANSFORMDATA transforms "y" data along ax1-ax2 2D plane
%
%Description: 
%    This function transforms "y" data along the 2D plane defined by ax1,ax2.
%    For example, let's say that ax1=1, ax2=size(E,1), rotAng = 27°. Then ax1 
%    is the major eigenvector, ax2 is the minor eigenvector. The major and 
%    minor eigenvectors become then 0° and 90° in the rotation plane. In 
%    consequence, 27° becomes a linear combination that falls somewhere in 
%    between. 
%    
%    For applying this transform we use a rotation matrix technique as to
%    rotate the major and minor vectors and forming a rotated basis matrix.
%    We then rotate by the original basis matrix, then apply the scaling, and 
%    revert by the new rotated basis.  More details are available in the  
%    thesis report Section 8.2.3, regarding Equation 8.24.
%
%Inputs: 
% 	 y - The residual data to process.
%    K - Scaling of two desirable dimensions. Null is [1,1]. 
%    rotAng - Rotation angle, in degrees. Null is 0. 
% 	 E - Eigenvectors matrix. This is the baseline of the transform.
%    ax1, ax2 - Indexes of eigenvectors to use as axes of the rotation plane.
%
%Outputs: 
% 	 y - The residual data after transforming. 
% 

% init.
doScale = (any(K~=1));
doRot = (rotAng~=0);
if (~doScale && ~doRot)
    return;
end

Erot = E;

% rotates requested axes
if (rotAng~=0)
    e1 = E(:,ax1);  % major axis
    e2 = E(:,ax2);  % minor axis
    ea1 = cosd(rotAng)*e1 + sind(rotAng)*e2;
    ea2 = -sind(rotAng)*e1 + cosd(rotAng)*e2;
    
    % create a new basis
    Erot(:,ax1)=ea1;
    Erot(:,ax2)=ea2;
end

% extra scaling by on 2 desirable dimensions
scl = eye(size(E));
if (doScale)
    scl(ax1,ax1) = K(1);
    scl(ax2,ax2) = K(2);
end
y = y * (E * scl * Erot');


function y = scalingData(y, E, ax1, ax2, scl)
%SCALINGDATA scales "y" data along two orthonormal dimensions.
%
%Description: 
%    This function scales "y" data along the two dimensions ax1,ax2. This 
%    is done by defining an eye matrix, except the two cells of the 
%    (ax1,ax2) pair. The data is then rotatins by the basis matrix, then 
%    the scaling is applied, and the rotation is revert by the same basis.
%
%Inputs: 
% 	 y - The residual data to process.
% 	 E - Eigenvectors matrix. This is the baseline of the transform.
%    ax1, ax2 - Indexes of eigenvectors to use as axes of the rotation plane.
%    scl - Scaling of two desirable dimensions. Null is [1,1]. 
%
%Outputs: 
% 	 y - The residual data after transforming. 
% 
%See also TRANSFORMDATA
%

if (all(scl==1))
    return;
end
SclMat = eye(size(E));
SclMat(ax1,ax1) = scl(1);
SclMat(ax2,ax2) = scl(2);
y = y * (E * SclMat * E');
