function rgbImg = cube2rgb(datacube)
%CUBE2RGB creates a representative RGB image for a multidimensional data. 
% 
%Description: 
% 	 This function applies PCA transform on a multi-layer image, and takes the 3-first
%    principal components. Since these components has the major eigenvalues and are 
%    orthogonal to each other, they can well represent the given data. The output
%    is normalized to [0,1] for visualization purposes. Note that this transform 
%    is also known as Karhunen-Loeve transform.
% 
%Inputs: 
% 	 datacube - A 3-D data-cube, of size [height, width, depth].
% 
%Outputs: 
% 	 rgbImg - An RGB image of the PCA visualization, of size [height, width, 3].
% 
%Example: 
%    figure; imshow(cube2rgb(data));  % creates a transformed data
% 
%See also PCA, PCA_CUBE
% 

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

data2_pca = pca_cube(datacube);
rgbImg(:,:,1) = imscale(data2_pca(:,:,1), minmax2(data2_pca(:,:,1)), [0,1]);
rgbImg(:,:,2) = imscale(data2_pca(:,:,2), minmax2(data2_pca(:,:,2)), [0,1]);
rgbImg(:,:,3) = imscale(data2_pca(:,:,3), minmax2(data2_pca(:,:,3)), [0,1]);
