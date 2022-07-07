function Create_Captions()
%CREATE_CAPTIONS creates figures and tables for the related paper.
%  
%Description: 
%    This script displays figures and tables by the order that they appear in the related paper:
%    Y. Furth, "Efficacy of Segmentation for Hyperspectral Target Detection". The script was run 
%    in MATLAB 9.5.0.944444 (R2018b), and image processing toolbox version 10.3 (R2018b). Basic 
%    concepts and conventions are detailed in "kernels\readme.txt". More background can be found  
%    in the related paper, and in the related thesis report, all available under
%    https://drive.google.com/drive/folders/13sYL6OAd45XWehQ0QEPVNDhR5WcMgGKC
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in gpl-3.0.txt).

% init.
Init_Env();
[d1, d2, d3] = Init_Data();

% Displays figures and tables by the order of the paper
Show_Figure1(d2.data, d2.seg, d2.t);
Show_Figure2(d1.data_orig, d2.data);
Show_Figure3(d1.seg, d2.seg);
Show_Figure4();
Show_Figure5();
Show_Figure11();
Show_Figure12();
Show_Figure13a();
Show_Figure13b(d1.data, d1.seg, d1.t);  
Show_Figure13c(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS, d2.E); 
Show_Figure14c(d1.data, d1.seg);  
Show_Figure15(d2.y, d2.seg, d2.phiG, d2.phiS, d2.E);  
Show_Figure16(d2.y, d2.seg, d2.phiG, d2.phiS, d2.E);  
Show_Figure17(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS);  
Show_Figure18(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS, d2.E, d3.wvlen);  
Disp_Table1a(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS);  
Disp_Table1b(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS, d2.E);  
Disp_Table1c(d2.data, d2.seg, d2.y, d2.phiG, d2.phiS);  
Disp_Table1d(d3.data, d3.seg, d3.t);  


function Init_Env()
%INIT_ENV adds required paths if needed
currDir = fileparts(mfilename('fullpath'));

fname = fullfile(currDir,'utils');
addPathIfNew(fname);

fname = fullfile(currDir,'kernels');
addPathIfNew(fname);


function addPathIfNew(fname)
%ADDPATHIFNEW adds path if not exist yet
pathCell = regexp(path, pathsep, 'split');
onPath = any(strcmpi(fname, pathCell));
if (~onPath)
    addpath(fname);
end
