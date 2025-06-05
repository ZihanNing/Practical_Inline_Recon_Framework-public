

function x = onesL(x, makeReal)

%ONESL creates an array with ones that has the same size and class as the input array. This avoids the repetitive use of ones(size(x),'like',x)
%
%   [X]=ONESL(X,{MAKEREAL})
%   * X is the input array.
%   * {MAKEREAL} to supress the use of complex array elements in case the input is complex.
%   ** X the output array with ones.
%
%   Yannick Brackenier 2022-01-30

if nargin < 2 || isempty(makeReal); makeReal = 0;end

N = size(x);
if makeReal
    x = ones( N , 'like', real(x));
else
    x = ones( N , 'like', x);   
end