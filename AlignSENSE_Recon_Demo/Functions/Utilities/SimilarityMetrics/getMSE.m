
function MSE = getMSE(x, xGT, ROI, useMagn)

%getMSE computes the Mean Squared Error (MSE) between an image and the ground truth (GT).
%   [MSE] = GETMSE(X,XGT,{ROI}) 
%   * X the array for which to compute the MSE.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** MSE the calculated MSE.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getMSE:: Ground Truth array expected.'); end
if nargin < 3 ; ROI = []; end
if nargin < 4 || isempty(useMagn); useMagn = 0; end
if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getMSE:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
if useMagn; x=abs(x);xGT=abs(xGT);end
errorImg = x - xGT;%Error of magnitude images

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%%% MSE
MSE = mean( abs(errorImg(:)).^2 );