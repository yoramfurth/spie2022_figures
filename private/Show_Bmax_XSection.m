function [h, tMax] = Show_Bmax_XSection(y, seg, yax1, yax2, a, th, constraintTMax, showKb, showBmax)
%SHOW_ELLIPSES_XSECTION plots Bmax over a cross-section as function of rotation angle.
% 
%Description: 
%    The axes of the cross-section is defined by 2 multi-dimensional non-parallel vectors (yax1, yax2).  
%    Over this plane, SHOW_ELLIPSES_XSECTION scans different angles (a), and plots the respective 
%    Bmax(a).
%
%Inputs: 
% 	 y - The residual data to process, generally formed by subtracting the local average from each pixel (x-m).
%    seg - Segmentation map. 
%    yax1, yax2 - The axes of the 2D cross-section.
%    a - Rotation angle of a target along planar axes.
%    th - The decision thresholds in terms of Pfa. 
%    constraintTMax - If true, grayouts angles out of the positive cone.
%    showKb - If true, shows also the graph of Kb(a).
%    showBmax - If true, plots annotates the location of the optimal Bmax(a). Grayed-out areas are ignored. 
%
%Outputs: 
% 	 h - Handles to the plotted objects.
% 	 tMax - The target that give the optimal Bmax(a). Grayed-out areas are ignored.
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Orthogonalizes and normalizes the axes (orthonormal projection)
yax1 = uvec(yax1);  % unit vector
yax2 = uvec(yax2);  % unit vector
yax2n = yax2 - yax1*(yax1'*yax2);
yax2n = uvec(yax2n);  

% Calculates Kb(a) values
aKb = a(1):0.1:a(end);
t_2d = [cosd(aKb); sind(aKb)];
[phiG, phiS] = Calc_PhiS_2d(y, seg, yax1, yax2);
Kb = SIM_Calc_Kb(phiG, phiS, t_2d);

% Calculates Bmax(a) values
Ntests = numel(a);
BMax = zeros(1,Ntests);
for n=1:Ntests
    fprintf('Ntest %d/%d\n', n, Ntests);
    t = sind(a(n)) * yax2n + cosd(a(n)) * yax1;
    t = uvec(t);    
    BMax(n) = Calc_Bmax(y, seg, 1, 0, t, th, 1); % fast estimation of Bmax 
end
fprintf('done!\n');

% Init. plots - defines constraint region
a_rng = [];
if (constraintTMax) 
    a_rng = Range_Positive_Combs(yax1, yax2n);  % constraints t to be positive 
end

% Finds the maximal Bmax(a) over all angles
tMax = [];
if (isempty(a_rng))
    nmax = argmax(BMax);
else
    ak = a+[0;360]; % a+360k
    kPos = find(any(ak>=a_rng(1) & ak<=a_rng(2)));
    nmax = kPos(argmax(BMax(kPos)));
end
if (nargout>=2)
    aMax = a(nmax);
    t = sind(aMax) * yax2n + cosd(aMax) * yax1;
    tMax = uvec(t);
end

% Plots Kb(a) graph on a left-axis 
legText = {};
h = [];
if (showKb)
    yyaxis left
    h(end+1) = plot(aKb, Kb);
    axis([a([1,end]),0,max(Kb(:))*2.5]);
    xlabel('\alpha [deg]'); ylabel('Kb'); grid on; hold on;
    legText = [legText,{'Kb'}];
    yyaxis right
end

% Plots Bmax(a) graph on a right-axis 
h(end+1) = plot(a, BMax);
axis([a([1,end]),0,max(BMax(:))*1.02]); 
xlabel('\alpha [deg]'); ylabel('Bmax'); grid on; hold on;
legText = [legText,{'Bmax'}];
ax=axis;

if (showBmax)
    h(end+1) = plot(a(nmax),BMax(nmax),'r-o');
    legTextTMax = 't_{max}';
    if (constraintTMax && ~isempty(a_rng))
        legTextTMax = 'max @t\geq0';
    end
    legText = [legText,{legTextTMax}];
end

% Grayouts the non positive regions
if (~isempty(a_rng))
    hf = GrayoutNegRegions(a_rng, ax);
    h = [h(:);hf(1)];
    legText = [legText,{'not  t\geq0'}];
end

% Plots dashed lines, representing main cross-section axes 
a0 = acosd(yax1'*yax2);
axesToPlot = [-180+[0,a0], 0,a0, 180+[0,a0]];
hpAxes = plot([1 1]'.*axesToPlot,ax(3:4)','--','Color',[1 1 1]*0.2);
h = [hpAxes(1); h(:)];
legText = [{'axes'},legText];
legend(h, legText,'location','SouthWest','FontSize',7);

% Arranges labels
hax = gca;
xticks = linspace(a(1), a(end), 9);
xticksLabel = split(strtrim(sprintf('%.3g° ', xticks)),' ');
set(hax,'XTickLabel',xticksLabel,'XTick',xticks);


function h = GrayoutNegRegions(a_rng, ax)
%GRAYOUTNEGREGIONS grayouts the non positive regions.

clrPosRng = [0.3,0.3,0.3];  % gray color
clrAlpha = 0.15;
if (a_rng(2)>180)
    a_rng = mod(a_rng+180, 360);
end
if (diff(a_rng)<0)
    a_neg = {a_rng([2,1])};  
else
    a_neg = {[-180 a_rng(1)], [a_rng(2) 180]};
end
fillNeg = @(a)fill(a([1 2 2 1]),ax([3 3 4 4]),clrPosRng,'EdgeColor','none','FaceAlpha',clrAlpha);
h = cellfun(fillNeg, a_neg);


function [phiGp, phiSp] = Calc_PhiS_2d(y, seg, yax1, yax2)
%CALC_PHIS_2D calculate the 2D covariances that appear over the cross-section plane.

% Estimates the local and global covariances (multi-dimensional)
[phiG, phiS] = Calc_phi(y, seg);

% Intersection with 2D plane
proj2d = [yax1, yax2]; % the 2D space (assuming perpendicular)
xs2dFun = @(phi) (proj2d' * (phi \ proj2d))^-1; % intersection of phi with 2d plane
phiGp = xs2dFun(phiG);
phiSp = cellfun(xs2dFun, phiS, 'UniformOutput', false);

