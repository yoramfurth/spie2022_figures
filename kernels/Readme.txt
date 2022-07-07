This folder contains functions related to the normalized segmented matched filter (NSMF), as introduced 
in the papers and the related thesis report. Below are detailed the main symbols and naming convention.


SOME BASICS:
------------
This library is designed for studying NSMF with specific context and assumptions. For example the target
is assumed to be smaller than a pixel, and to satisfy the additive model, that is, the yield for the 
"corrupted" pixel (x) is given by just appending the target's spectrum (t) with some small power (p),
that is, x:=x+pt. However, the functions were written modular and context-less enough to be easily reuse. 

The NSMF is a type of matched-filter that deals with segmentation assuming target additive model. Its 
definition is: 
        NSMF = t' * phi^-1 * (x-m) / sqrt(t' * phi^-1 * t)   ,
while "m" is the local average of each pixel, and "phi" is the covariance matrix of each segment (phiS),
alternately. The global version, sometimes called NGMF, is calculated by substituting phiG instead. 
Here we normally substitute y=x-m, which is solved preliminarily as a separate problem. The covariance 
"phi" therefore is estimated estimated directly from "y" by applying y*y', that is, correlation of 
any pair of channels. 

More details are introduced in the related paper, and in the related thesis report ("Efficacy of
Segmentation for Hyperspectral Target Detection", Yoram Furth 30/9/2021), all available under
https://drive.google.com/drive/folders/13sYL6OAd45XWehQ0QEPVNDhR5WcMgGKC


NAMING CONVENTION:
------------------
The function names are built of: the context, the type of operation, the objective. For example: 
"SIM_Calc_mu" means a whole calculation of E(NSMF) on synthetic data. This works as follows.

The context - There are 2 types of data to be operated: synthetic for simulations, and real data-cubes 
for diverse purposes. These leads to significantly different calculations, so different functions. 
For differing them, we added a "SIM_" suffix to scripts dedicated especially for the former case.

Type of operation - Dedicated keywords are used for hinting the type or operation. The major of them
are as follows. "Calc" refers to the whole calculation of a term, "Eval" refers to a partial calculation, 
that is, given another term\s, already calculated. 

The objective - Generally on of the symbols here below. 


SYMBOLS NAMING:
---------------
x, data - Hyperspectral data 
m - Estimated background
y - Estimated residual noise, generally (x-m)
s - Indexes of spectral clusters 
phi, phiG, phiS - Covariance matrix / global / local per-segment
t -  A vector representing the spectrum of the target of interest (often refers to its spectral direction only)
p - Portion of a target within a pixel (often refers to the total target's power)

z - Detector score axis (used also as a discriminant axis)
fN / fT - PDF of the scores of NSMF without / with target
cN / cT - CDF of the scores of NSMF without / with target
eta / etaG / etaL - Decision threshold on z-axis / in global domain / in local domain
mu / muG / mu1,mu2 - Score expectation / in global domain / per segment (usually "with target")
q / qG / qL - Score expectation when p=1, which is equivalent to SNR in the NSMF case / evaluated globally / locally (an array)
pG / pL - Special p-anchors defined by eta/q, evaluated globally / locally (an array) 

r1,r2 = Axes of ellipse (corresponds to "a,b" in the documents)
AR - The ellipse aspect-ration (r1/r2), often coming from an ellipsoid major and minor axes
K - Scalar multiplier that represents a data rescaling amount
theta - Opening angle between a pair of ellipsoid's majors [degrees]
a - Rotation angle of a target along planar axes (corresponds to "alpha" in the documents)

Pfa / Pd - Probability of false alarm / detection above "eta" threshold
th - Decision threshold in terms of Pfa
A / AG / AL - Detection algorithm evaluator as a function of "th" / in global / in local domain
B - Segmentation benefit evaluator as a function of a specific detector's factors

pMax, Bmax - The optimum of B(p) function
tmax - The optimum of B(t) function
Ksplit / Kshift / Kb - Directional non-stationarity impact of type split / shift / both


ABBREVIATIONS:
--------------
MF / GMD / SMF - Matched filter / global matched filter / segmented matched filter
NGMF / NSMF - Normalized global / segmented matched filter
PCA - Principal components analysis, related to Karhunen–Loeve transform
ROC - Receiver operating characteristic curve
SIC - Special in-house real data cube, used for simulations
AUC - Area under ROC curve, normally in a range of thresholds
SNR - Signal to noise ratio


COPYRIGHT:
----------
Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
Dept. Electrical & Computer Engineering, BGU Israel.
This code is published under GNU GPLv3 license (see license in gpl-3.0.txt).
