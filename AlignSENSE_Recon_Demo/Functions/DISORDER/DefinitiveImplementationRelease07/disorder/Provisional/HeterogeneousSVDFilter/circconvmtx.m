function y=circconvmtx(x,N,c)

%CIRCCONVMTX  creates the circular convolution matrix of size N for x (as 
%a row vector)
%   Y=CIRCCONVMTX(X,N,C)
%   * X is the convolution filter
%   * {N} is the space of the signals to convolve. It defaults to the
%   length of X
%   * {C} serves to center the convolution. It defaults to 0
%   ** Y is the convolution kernel
%

if nargin<3 || isempty(c);c=0;end

x=x(:);
M=size(x,1);
if nargin<2;N=M;end
assert(N>=M,'Size of space assumed bigger than size of kernel');
x(end+1:N,1)=0;
y=toeplitz(x,[x(1);flip(x(2:N))]);%Not checked whether conjugation is needed here
if c;y=circshift(y,[0 floor((M-1)/2)]);end
