function Show_Figure2(SIC_data, COOKE_data)
%SHOW_FIGURE2 shows Figure 2 - datasets used for this study.
%
%Inputs: 
% 	 SIC_data - SIC data-cube. 
% 	 COOKE_data - Cooke's data-cube. 
%

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% Jet colormap transformed from jet to gray (more clear than a simple gray display)
cmapGrayJet = imscale(rgb2gray(reshape(jet(64),[1,64,3])), [5/255, 240/255], [0 1])' * [1,1,1];

% Figure 2(a) - shows band 15 from SIC data
nBand = 15;
figure; imagesc(SIC_data(:,:,nBand))
set(gcf,'position',[831, 563, 551, 408], 'colormap', cmapGrayJet)

% Figure 2(b) - shows band 50 from Cooke's data
nBand = 50;
figure; imagesc(COOKE_data(:,:,nBand))
set(gcf,'position',[651, 504, 1269, 427], 'colormap', cmapGrayJet)
