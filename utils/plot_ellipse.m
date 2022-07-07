function hPlot = plot_ellipse(r1, r2, x0, y0, phi, varargin)
%PLOT_ELLIPSE draws an ellipse according to its components. 
%
%Description: 
%    This function plots on the current axes an ellipse that represents the requested
%    parameters, including its radius its centers, and the angle of rotation. 
%
%Inputs: 
%    r1,r2 - Value of the major and minor axis, respectively.
%    x0,y0 - Abscissa and Ordinate of the center point of the ellipse, respectively.
%    phi - Angle [radians] between x-axis and the ellipse major axis. 
%    varargin - Parameters to pass to plot() function, e.g. line style.
%    
%Outputs: 
%    hPlot - A handle to line object.
%
%Example: 
%    % plots an ellipse with size of 5,3, around 1,-2, rotated by pi/4
%    h = plot_ellipse(5,3,1,-2,pi/4,'r--');  % linestyle if 'r--'
%    set(h,'LineWidth',2);  
%
%See also PLOT
%

%    Copyright 25/3/2003 Lei Wang (WangLeiBox@hotmail.com), as "ellipsedraw.m"
%    Dept. Mechanical & Aerospace Engineering, NC State University.
%    https://www.mathworks.com/matlabcentral/fileexchange/3224-ellipsedraw1-0
%    
%    Updated by Yoram Furth, 24/06/2017, reorganized and renamed to "plot_ellipse.m".
%    Dept. Electrical & Computer Engineering, BGU Israel.

if (nargin<5)
    phi = 0;
end
if (nargin<4)
    x0 = 0; y0 = 0;
end

theta = (0.0001: 0.01 :2*pi);

% Parametric equation of the ellipse
x = r1 * cos(theta);
y = r2 * sin(theta);

% Coordinate transform
X = cos(phi) * x - sin(phi) * y;
Y = sin(phi) * x + cos(phi) * y;
X = X + x0;
Y = Y + y0;

% Plot the ellipse
hPlot = plot(X, Y, varargin{:});
