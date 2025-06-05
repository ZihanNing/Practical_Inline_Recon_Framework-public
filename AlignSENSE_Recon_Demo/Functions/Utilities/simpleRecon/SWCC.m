

function x = SWCC (y, S)

%SWCC computes the Sensitivity Weighted Coil Combination (SWCC) of multi-channel
%data in the image domain.
%   [X] = SWCC(Y,S) 
%   * Y the array containing multi-channel image data.
%   * S the array containing coil sensitivities of every channel.
%   ** X the SWCC image.
%
%   Yannick Brackenier 2022-01-30
    
NY = size(y);NS = size(S);
dimChannel = numDims(S);

assert(isequal(NY(1:3),NS(1:3)),'SWCC: Array size of y and S need to match.');

x = bsxfun(@times,...
    multDimSum(  bsxfun(@times,conj(S),y),  dimChannel),...
    1./(normm(S,[],dimChannel)+1e-6));

end