function x=signz(x)

%SIGNZ   Calculates the sign function under the convention that zero values
%return value 1 instead of 0
%   X=SIGNZ(X)
%   * X is the input array
%   ** X is the output array
%

x=sign(x);
x(x==0)=1;
