

function [SNR] = getSNR(x,xGT,ROI,useMagn)

%GETSNR computes the signal to noise ratio (SNR) between an image and the ground truth (GT).
%   [SNR] = GETSNR(X,XGT,{ROI}) 
%   * X the array for which to compute the SNR.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** SNR the calculated SNR in dB.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getSNR:: Ground Truth array expected.'); end
if nargin < 3 ; ROI = []; end
if nargin < 4 || isempty(useMagn); useMagn = 0; end
if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getSNR:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
if useMagn; x=abs(x);xGT=abs(xGT);end
errorImg = x - xGT;   

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%%% Compute supplementary info
sm2 = multDimMea(normm(x)); %Means Squared Signal
%var = std([ real(errorImg(:)) ; imag(errorImg(:)) ] )^2; %Variance of real and imaginary component (stacked) - DEPRECATED
var =  multDimMea(normm(errorImg)); %Means Squared Signal of error

%%% SNR
SNR = gather( 10*log10(sm2 / var) ) ;
