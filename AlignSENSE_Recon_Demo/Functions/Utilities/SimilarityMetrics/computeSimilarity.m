
function [SI] = computeSimilarity(x,xGT,type,ROI,useMagn)

%COMPUTESIMILARITY computes a similarity metric between an image and the corresponding ground truth (GT).
%   [SI] = COMPUTESIMILARITY(X,XGT,{TYPE},{ROI}) 
%   * X the array for which to compute the metric.
%   * XGT the ground truth array.
%   * {TYPE} the type of similarity metric to compute.
%   * {ROI} the region of interest (ROI) for which to compute the similarity metric. It must be a binary mask. 
%   * {USEMAGN} whether to compute the metric on the magnitude of the images (in case phase is unreliable). 
%   ** SI the calculated similarity metric.
%        for type 'SNR' --> In dB
%        for type 'pSNR' --> In dB
%        for type 'AP' --> In dB
%        for type 'MI' --> In a.u.
%        for type 'SSIM' --> In a.u.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(xGT); error('computeSimilarity:: Ground Truth array expected.'); end
if nargin < 3 || isempty(type); type = {'SNR'}; end
if nargin < 4 || isempty(ROI); ROI = []; end
if nargin < 5 || isempty(useMagn); useMagn = 0; end

if ~ isempty(ROI); assert(all(ismember(ROI(:), [0 1])),'computeSimilarity:: Must be binary mask. Weighted similarity to be implemented!');end
    

SI = zeros( [1,length(type)] );

for i=1:length(type)
    
    if strcmp( type{i}, 'SNR');         SI(i) = getSNR(x, xGT, ROI, useMagn);
    elseif strcmp( type{i}, 'pSNR');    SI(i) = getpSNR(x, xGT, ROI, useMagn);
    elseif strcmp( type{i}, 'SSI');    SI(i) = getSSI(x, xGT, ROI, useMagn);
    elseif strcmp( type{i}, 'MI');      SI(i) = getMI(x, xGT, ROI);%MI by default computed on magnitude images
    elseif strcmp( type{i}, 'AP');      SI(i) = getAP(x, xGT, ROI, useMagn);    
    elseif strcmp( type{i}, 'MSE');     SI(i) = getMSE(x, xGT, ROI, useMagn);
    end
    
end
