
function x = RSOS (y,dim)

%RSOS computes the Root Sum Of Squares (RSOS) combination of multi-channel
%data in the image domain.
%   [X] = RSOS(Y) 
%   * Y the array containing multi-channel image data.
%   * {DIM} the dimension along which to combine channels.
%   ** X the RSOS image.
%
%   Yannick Brackenier 2022-01-30
    
if nargin<2 || isempty(dim);dim=4;end

x = sqrt(normm(y,[],dim));

end