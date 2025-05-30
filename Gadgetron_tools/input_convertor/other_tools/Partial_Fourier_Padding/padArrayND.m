
function x = padArrayND (x,padDim,di,padValue, padType)

%PADARRAYND pads an N-dimensional array.
%   [X]=PADARRAYND(X, PADDIM {DI},{PADVALUE},{PADTYPE})
%   * X array to pad.
%   * PADDIM the number of elements to pad to each dimension.
%   * {DI} the direction of padding. 0 means extracting the original array from the padded array.
%   * {PADVALUE} the value to pad (Defaults to the first value of the array).
%   * {PADTYPE} to type of padding ('both','pre','post').
%   ** X is the padded array.
%
%   Yannick Brackenier 2022-01-30

if nargin <2 || isempty(padDim); return; end
if nargin <3 || isempty(di); di = 1;end
if nargin <4 || isempty(padValue); padValue = gather(x(1));end %Background value. padarray does not support gpu array for padValue.
if nargin <5 || isempty(padType); padType = 'both';end%Can be 'pre' or 'post'

N =size(x);
NDims = numDims(x);

if length(padDim) < NDims; padDim(end+1:NDims)=0;end
    
if di %Forward padding
    x = padarray(x,double(padDim), padValue, padType);

else %ROI extraction   
    ROIDyn=cell(1,NDims);
    for n=1:length(ROIDyn);ROIDyn{n}=(padDim(n)+1):(N(n)-padDim(n));end

    x = dynInd(x,ROIDyn,1:NDims) ;     
end

