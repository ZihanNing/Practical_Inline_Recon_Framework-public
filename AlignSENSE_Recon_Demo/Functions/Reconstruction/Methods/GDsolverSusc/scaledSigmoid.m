

function y = scaledSigmoid(x, a)

%x input
%a scaling in variable direction

if nargin<2 || isempty(a); a =1;end


y = 1 ./ (1 + exp(-a* x) );