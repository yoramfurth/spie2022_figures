function Show_Figure3(SIC_seg, COOKE_seg)
%SHOW_FIGURE3 shows Figure 3 - Segmentation maps used for the datasets.
%
%Inputs: 
% 	 SIC_seg - SIC data-cube's segmentation (2 segments). 
% 	 COOKE_seg - Cooke's data-cube's segmentation (5 segments).
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Figure 3(a) - shows SIC data-cube's segmentation
figure,imagesc(SIC_seg)
set(gcf,'position',[797,695,304,187],'colormap',[0.5,0,0; 0,0,0.5625])
colorbar('Ticks',[1.25,1.75],'TickLabels',{'1','2'})

% Figure 3(b) - shows Cooke's cube's segmentation
figure,imagesc(COOKE_seg)
set(gcf,'position',[1071,87,800,280],'colormap',parula(5))
colorbar;
