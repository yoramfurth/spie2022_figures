function hPlot = plotline(k,hAx,lineStyle,lineWidth,color)
%PLOTLINE plots vertical lines on the given axes.
%
%Description: 
%    This function draws vertical line with automatic height at the given
%    points along the horizontal axis. 
%
%Inputs: 
%    k - index on horizontal axis. Supports array of multiple values.
%    hAx - handle to the given axes.
%    lineStyle, lineWidth, color - parameters for plot() built-in function.
%
%Outputs: 
%    hPlot - a handle to line object.
%
%Examples: 
%     y1 = 1;
%     y2 = 3; 
%     PLOTLINE([y1 y2],gca,2,'r') % plots 2 vertical lines on the current axis.
%
%See also: plot
%

%    Copyright 2019-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
marginFactor = 1.1;
ax = axis(hAx);  % remind axis

% extracts graphs from figure
hLine = findobj(hAx,'Type','line');
lineLims = min(get(hAx,'YLim'));

ft = nan(numel(hLine),2);

for kh = 1:numel(hLine)
    xLine = get(hLine(kh),'XData');
    yLine = get(hLine(kh),'YData');
    if ((numel(xLine))<3)
        continue;
    end
    
    yMax = max(yLine(:));
    ft(kh,:) = [yMax, yMax];
end

% ignore nan points
ftMax = max(ft,[],1);
if (all(isnan(ftMax)))
    hPlot = 0;
    return;
end
ftMax(isnan(ftMax)) = mean(ftMax(~isnan(ftMax)));

% determines lines to plot
xx = []; yy = [];
for kt=1:numel(k)
    lineLimsMax = max(ftMax(:)); %ftMax(kt);
    lineLims(2) = lineLims(1) + (lineLimsMax - lineLims(1)) * marginFactor;
    if (kt>1)
        hold on;
    end
    xx = [xx; k([kt kt])]; %#ok
    yy = [yy; lineLims]; %#ok
end

% plot the line
hPlot = plot(xx', yy', lineStyle, 'linewidth', lineWidth, 'Color', color);

% retrieves the original axis
axis(hAx,ax);

