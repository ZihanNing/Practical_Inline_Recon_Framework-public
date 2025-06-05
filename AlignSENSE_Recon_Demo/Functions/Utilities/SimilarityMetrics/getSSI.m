function [SSI] = getSSI(x,xGT, ROI, useMagn, multiScale)

%GETSSI computes the Structural Similarity Index (SSI) between an image and the ground truth (GT). The SSI is described in Usman, M., Latif, S., Asim, M., Lee, B. D., & Qadir, J.(2020). Retrospective Motion Correction in Multishot MRI using Generative Adversarial Network. Scientific Reports, 10(1). https://doi.org/10.1038/s41598-020-61705-9.
%   [SSI] = GETSSI(X,XGT,{ROI},{MULTISCALE}) 
%   * X the array for which to compute the SSI.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   * {MULTISCALE} a flag to use the multi-scale version of the SSI implemented in matlab. (Defaults to 0).
%   ** SSI the calculated SSI.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getSSI:: Ground Truth array expected.'); end
if nargin < 3 ; ROI = []; end
if nargin < 4 || isempty(useMagn);  useMagn = 0; end
if nargin < 5 ; multiScale = 0; end

if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getSSI:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
if useMagn; x=abs(x);xGT=abs(xGT);end
errorImg = x - xGT;   

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%Gather as SSI not defined for GPU arrays
x= gather(x);
xGT=gather(xGT);

if multiScale
    SSI = multissim3(x,xGT);
else
    SSI = ssim(x,xGT);
end