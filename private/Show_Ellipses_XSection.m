function h = Show_Ellipses_XSection(phiG, phiS, yax1, yax2, t, constraintTMax)
%SHOW_ELLIPSES_XSECTION plots ellipsoids cross-section.
%
%Description: 
%    A cross-section of an ellipsoid gives a 2D ellipse. This function plots these ellipses
%    for the ellipsoids that represent covariance matrixes. The cross-section plane is 
%    represented by two non-parallel vectors (yax1, yax2). 
%
%Inputs: 
%    phiG - The global covariance matrix that corresponds to all the local {phiS} together.
%    phiS - Cell array of local covariance matrixes (i.e. per segment).
%    yax1, yax2 - The axes of the 2D cross-section.
%    t - A vector that represents a target's spectrum (optional).
%    constraintTMax - If true, grayouts angles out of the positive cone.
% 
%Outputs: 
% 	 h - Handles to the plotted objects.
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Orthogonalizes and normalizes the axes (orthonormal projection)
yax1 = uvec(yax1);
[yax2, yvec2] = deal(uvec(yax2));
corr12 = yax1'*yax2;
if (abs(corr12)>1e-4)  
    yax2 = yax2 - yax1*corr12;
    yax2 = uvec(yax2);
end

% Projects everything on the 2D plane
proj2d = [yax1, yax2]; % the 2D space
xs2dFun = @(phi) (proj2d' * (phi \ proj2d))^-1; % intersection of phi with 2d plane
phiGp = xs2dFun(phiG);
phiSp = cellfun(xs2dFun, phiS, 'UniformOutput', false);

projVec2d = (proj2d'*proj2d) \ proj2d';
tmaxp = projVec2d * t;
tmaxp = tmaxp * (norm(tmaxp)/norm(t));
pMaxEst = sqrt(2) / max(cellfun(@(phi) Eval_q(tmaxp, phi), phiSp));

% Init. plots - defines constraint region
a_rng = [];
if (constraintTMax)
    a_rng = Range_Positive_Combs(yax1, yax2);  % constraints t to be positive 
end

% Plots ellipses
h = Show_Ellipses(phiSp, phiGp, tmaxp*pMaxEst);
vec2legendcell = @(v) split(strtrim(sprintf('\\Phi_%d ', v)),' ');
legText = [vec2legendcell(1:numel(phiS))', {'\Phi_G'}];
legTextTMax = 't_{max}';
if (constraintTMax && ~isempty(a_rng))
    legTextTMax = 'max @t\geq0';
end
legText = [legText, {legTextTMax}];

% Grayouts the non positive region
if (~isempty(a_rng))
    clrPosRng = [1 1 1]*0.3;     % gray
    clrAlpha = 0.15;
    ax = axis;
    ta_rng = cosd(a_rng) .* [1;0] + sind(a_rng) .* [0;1];
    polyAx = polyshape(ax([1 1 2 2]),ax([3 4 4 3]));
    polySlice = polyshape([ta_rng' * 4 * max(ax(:)); 0 0]);
    polyNeg = subtract(polyAx, polySlice);
    h(end+1) = plot(polyNeg,'FaceColor',clrPosRng,'EdgeColor','none','FaceAlpha',clrAlpha);
    legText = [legText, {'not  t\geq0'}];
end

% Plot axes bars
yaxesToPlot = [-yax1, -yvec2, yax1, yvec2];  % prepare axes plotting
axFactor = max(axis)*2;
yax2p = projVec2d * yaxesToPlot * axFactor; % the 2D space
orig = zeros([1,size(yax2p,2)]);
hp = plot([orig;yax2p(1,:)],[orig;yax2p(2,:)],'Color',[1 1 1]*0.4);
uistack(hp, 'bottom');

% Adds legend
h = [hp(1);h];
legText = [{'axes'},legText];
legend(h, legText, 'location', 'NorthWest','FontSize',7)


function h = Show_Ellipses(phiS, phiG, t)
%SHOW_ELLIPSES draws 2D ellipses of covariance matrixes and a vector of a target.

% Init.
N = numel(phiS);
h = zeros(N+1,1);
axis equal
hold on;

% Plots everything: phiS{}, phiG, t
for n=1:N
    h(n) = cov_ellipse_draw(phiS{n});
end
h(N+1) = cov_ellipse_draw(phiG,'k--','LineWidth',2);
h(N+2) = SIM_Show_t(t*1.05, 'Color', [1,0,0]); % additional factor for pretty plot

% Arranges the figure
ax = [-1 1 -1 1]*max(axis)*1.1;
axis(ax);
grid on;

vec2legendcell = @(v) split(strtrim(sprintf('\\Phi_%d ', v)),' ');
legText = [vec2legendcell(1:N)',{'\Phi_G','t'}];
legend(legText,'location','NorthWest')


function h = cov_ellipse_draw(Sig,varargin)
%COV_ELLIPSE_DRAW draws a 2D ellipse that corresponds to a covariance matrix.

% The center
Mu = [0 0];
x0 = Mu(1);
y0 = Mu(2);

% Radius and rotation angle
[V,D,~] = svd(Sig);
phi = atan2(V(2,1),V(1,1));
r1 = sqrt(D(1,1));
r2 = sqrt(D(2,2));

% Plot
h = plot_ellipse(r1, r2, x0, y0, phi, varargin{:});


function h = SIM_Show_t(t,varargin)
%SIM_SHOW_T draws a 2D vector that represents a target's spectrum (t).

% Sets default color as red
if (~any(strcmpi(varargin,'color')))
    varargin = [varargin(:)',{'Color', [1,0,0]}];
end

% Plot the vector
h = plot([0 t(1)], [0 t(2)], 'LineWidth', 1.5, varargin{:});
hold on; plot(t(1), t(2), 'o', 'LineWidth', 1.5, varargin{:});
