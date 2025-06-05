
function MI = getMI(x, xGT, ROI)

%GETMI computes the Mutual Information (MI) between an image and the ground truth (GT).
%   [MI] = GETMI(X,XGT,{ROI}) 
%   * X the array for which to compute the MI.
%   * XGT the ground truth array.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** MI the calculated MI.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('getMI:: Ground Truth array expected.'); end
if nargin < 3 ; ROI = []; end
if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'getMI:: Must be binary mask. Weighted SNR to be implemented!');end

%%% Compute error
errorImg = abs(x) - abs(xGT);%Error of magnitude images

%%% Restrict to ROI
if ~isempty(ROI)
    errorImg(ROI==0) = []; 
    x(ROI==0)=[]; 
    xGT(ROI==0)=[];
end

%Gather as MI not defined for GPU arrays
x= gather(x);
xGT=gather(xGT);

%%% MI
MI = mi(x,xGT);
