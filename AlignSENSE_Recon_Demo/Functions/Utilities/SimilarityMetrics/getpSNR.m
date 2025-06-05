
function [pSNR] = getpSNR(x,xGT,ROI,useMagn)

%GETPSNR computes the peak SNR (pSNR) between an image and the ground truth (GT).Described in Usman, M., Latif, S., Asim, M., Lee, B. D., & Qadir, J.(2020). Retrospective Motion Correction in Multishot MRI using Generative Adversarial Network. Scientific Reports, 10(1). https://doi.org/10.1038/s41598-020-61705-9.
%   [PSNR] = GETPSNR(X,XGT,{ROI}) 
%   * X the array for which to compute the PSNR.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** PSNR the calculated PSNR in dB.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getpSNR:: Ground Truth array expected.'); end
if nargin < 3 || isempty(ROI); ROI = []; end
if nargin < 4 || isempty(useMagn); useMagn = 0; end
if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getpSNR:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
if useMagn; x=abs(x);xGT=abs(xGT);end
errorImg = x - xGT;   

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%%% SNR
pSNR= 20*log10( max(abs(errorImg(:))) /... 
                sqrt(mean(abs(errorImg(:)).^2)));
