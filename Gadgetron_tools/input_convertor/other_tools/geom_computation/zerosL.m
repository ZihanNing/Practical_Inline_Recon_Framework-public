
function x = zerosL(x, makeReal)

%ZEROSL creates an array with zeros that has the same size and class as the input array. This avoids the repetitive use of zeros(size(x),'like',x)
%
%   [X]=ZEROSL(X,{MAKEREAL})
%   * X is the input array.
%   * {MAKEREAL} to supress the use of complex array elements in case the input is complex.
%   ** X the output array with zeros.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(makeReal); makeReal = 0;end

N = size(x);
if makeReal
    x = zeros( N , 'like', real(x));
else
    x = zeros( N , 'like', x);   
end