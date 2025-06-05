
function AP = getAP(x, xGT, ROI, useMagn)

%GETAP computes the artefact power (AP) between an image and the expected ground truth (GT). The AP is based on Omer, H. & Dickinson, R. A graphical generalized implementation of sense reconstruction using MATLAB. Concepts in Magnetic Resonance Part A 36(3), 178–186 (2010).
%Also described in Usman, M., Latif, S., Asim, M., Lee, B. D., & Qadir, J.(2020). Retrospective Motion Correction in Multishot MRI using Generative Adversarial Network. Scientific Reports, 10(1). https://doi.org/10.1038/s41598-020-61705-9.
%   [AP] = GETAP(X,XGT,{ROI}) 
%   * X the array for which to compute the AP.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** AP the calculated AP in dB.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getAP:: Ground Truth array expected.'); end
if nargin < 3 ; ROI = []; end
if nargin < 4 || isempty(useMagn); useMagn = 0; end
if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getAP:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
if useMagn; x=abs(x);xGT=abs(xGT);end
errorImg = abs(x) - abs(xGT);%Error of magnitude images

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%%% Aretefact Power
AP = sum( errorImg(:).^2)/ sum(xGT(:).^2);